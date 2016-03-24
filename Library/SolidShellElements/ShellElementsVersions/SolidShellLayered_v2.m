classdef SolidShellLayered_v2 < handle
    
    %Version 2 of layered solid shell element.
    %One EAS-M-matrix for ALL layers
    %One combinded (higher order) interpolation for the stresses
    
    
    properties
        
        %Layered integration rule
        lir;% = IntegrationRule;
        
        %Interolationrules
        dispInterp;
        stressInterp;
        stressInterpOrder; %interpolation order for each stresss
        
        submatrices;
        
        %Interpolation matrix for strains
        Mhat = -1;
        
        %Number of dofs
        nDispDofs = -1;
        nStrainDofs = -1;
        nStressDofs = -1;
        nLamStrainDofs = -1;
        
        %Element coords
        ez; ey; ex;
        
        %Number of laminas
        nLam;
        
        %The deterimatant of each layer (always concstant for each element)
        lay_detJ;
        
        %Dont know if alpha and beta should be stored in the element, by I
        %do it here temporararly
        beta;
        alpha; %layer-wise
    end
    
    methods
        function obj = SolidShellLayered_v2(ngpx, ngpy, ngpz, ex, ey, ez, stressInterp, Mhat, elprop)
            
            %Number of layers
            obj.nLam = elprop.nLam;
            
            %Create integration rules
            obj.lir = LayeredIntegrationRule(elprop.nLam, ngpx, ngpy, ngpz);
            
            %Create Strain interpolation matrices
            obj.Mhat = Mhat;
            
            %Define interpolation for each stress component
            for ii=1:length(stressInterp)
                val = stressInterp(ii);
                if(val == 2)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z2;
                elseif(val == 3)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z3;
                elseif(val == 4)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z4;
                elseif(val == 7)
                    obj.stressInterp{ii} = InterpolatorX2Y2Z7;
                else
                    error('Interpolator not implemented yet.');
                end
            end
            obj.stressInterpOrder = stressInterp;
            
            %Define number of dofs
            obj.nStrainDofs = size(Mhat(0,0,0),2);
            obj.nStressDofs = sum(stressInterp*2*2);
            obj.nDispDofs = 3*2*2*2; %*nnoz??
            
            %Save element coords
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
            %Interpolation for displacement always the same
            obj.dispInterp = InterpolatorX2Y2Z2;
            
        end
        
        function [R, J] = computeRandJ(obj,a, eq, ftrac, elprop)
            
            %Init
            V  = 0;
            Ae = zeros(obj.nStrainDofs, obj.nDispDofs);
            Ce = zeros(obj.nStrainDofs, obj.nStrainDofs);
            Le = zeros(obj.nStressDofs, obj.nStrainDofs);
            Qe = zeros(obj.nStressDofs, obj.nStressDofs);
            Ke = zeros(obj.nDispDofs, obj.nStressDofs);
            Re = zeros(obj.nStressDofs, obj.nDispDofs);
            fe = zeros(obj.nDispDofs,1);

            %Calculate alpha values
            for ilay = 1:obj.nLam
                layD = elprop.Dmatrices(:,:,ilay);             

                for gp = obj.lir.irs(ilay).gps %obj.lir.irs(ilay).gps
                    
                    %Gauss coordinates in local-layer-system and
                    %local-element-system (only need to map the z-coord
                    %since the x and y-coords should be the same
                    lcoords = gp.local_coords;
                    ecoords = lcoords;
                    ecoords(3) = obj.lir.getElementGaussCoordinate(...
                        lcoords(3), elprop.int_coordsL(ilay:ilay+1));
                    
                    %N-vector, eval in element-system
                    Nvec = obj.dispInterp.eval_N(ecoords);
                    
                    %Derivatives of n-vector
                    [dNdx, detJEl] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
                    
                    % Get detJ for the layer
                    layZ = [[1,1,1,1]*elprop.int_coordsG(ilay), [1,1,1,1]*elprop.int_coordsG(ilay+1)];
                    [~, detJ] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, layZ);
                    
                    % Get N and B matrix used in FEM
                    [N, B] = solid8NandBmatrix(Nvec, dNdx);
                    
                    % Enanced part, eval in element-system
                    Mi = obj.Mhat(ecoords(1), ecoords(2), ecoords(3));
                    
                    %Stress part, eval in element-system
                    P = [];
                    % Varför inte är inte alla 6 spänningskomponenter med? /JB
                    for is = 1:length(obj.stressInterp)
                        NvecSigma = obj.stressInterp{is}.eval_N(ecoords);
                        % TODO: Mycket långsamt med blkdiag! /JB
                        P = blkdiag(P,NvecSigma);
                    end
                    
                    % Integration
                    dV = detJ * gp.weight;
                    
                    V = V + 1*dV;
                    Ae = Ae + Mi'*layD*B   * dV;
                    Ce = Ce + Mi'*layD*Mi  * dV;
                    Le = Le + P' *layD*Mi  * dV;
                    Re = Re + P' *layD*B   * dV;
                    Qe = Qe + P'*P         * dV;
                    Ke = Ke + B'*P         * dV;
                    
                    fe = fe + N'*eq   * dV;
                    
                end
            end

            %add the traction vector to the total force vector.
            fext = fe + ftrac;
            alpha = -Ce\(Ae*a);
            
            %Want to solve: ...Qe\(Re + Lsum);..., but we have some bc:s
            %Bc on sigma
            [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
            lockedSigma = [TxybDofs, TxytDofs]';
            sigmaBc = [lockedSigma, lockedSigma*0];
            lockedSigma = [TxybDofs, TxytDofs]';
            allSigma = 1:obj.nStressDofs;
            freeSigma = setdiff(allSigma, lockedSigma);

            %Solve for beta and derivate of beta wrt a
            %Caluclate righthandside of equations
            dfbeta = Re - Le*inv(Ce)*Ae;
            fbeta  = dfbeta*a;
            
            %Init beta and dbeta/da
            beta    = zeros(obj.nStressDofs,1);
            dbetada = zeros(obj.nStressDofs,obj.nDispDofs);

            %Solve (dont forget the bc:s :-) )
            beta   (freeSigma  ) = Qe(freeSigma,freeSigma)\(fbeta(freeSigma) - Qe(freeSigma,lockedSigma)*sigmaBc(:,2));
            dbetada(freeSigma,:) = Qe(freeSigma,freeSigma)\(dfbeta(freeSigma,:));

            R = Ke*beta - fext;
            J = Ke*dbetada;

            %Save internal variables for element
            obj.beta = beta;
            obj.alpha = alpha;

        end
        
        function dofs = getStressComponentDofs(obj, component, nodes)
            
            ss = obj.stressInterpOrder;
            
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i)*2*2 ;
                startFrom = startFrom + temp;
            end
            
            dofs = startFrom + nodes;
            
        end
        
        function [bottomDofs, topDofs] = getStressComponentTopAndBottonDofs(obj, component)
            
            ss = obj.stressInterpOrder;
            
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i)*2*2 ;
                startFrom = startFrom + temp;
            end
            
            temp = 1:obj.stressInterpOrder(component)*2*2;
            
            dofs = startFrom + [temp(1:4), temp((end-3):end) ];
            bottomDofs = dofs(1:4);
            topDofs    = dofs(5:8);
        end
        
        function P = getStressPmatrix(obj, lcoords)
            
            P = [];
            for is = 1:length(obj.stressInterp)
                NvecSigma = obj.stressInterp{is}.eval_N(lcoords);
                P = blkdiag(P,NvecSigma);
            end
            
        end
    end
    
end

%
%             allSigma = 1:obj.nSigmaDofs;
%             freeSigma   = setDiff(allSigma, lockedSigma)
%
%             %Bc vector
%             sigmaBc = [lockedSigma, lockedSigma*0];
%
%             %Get new "stress-matrices"
%             RL = (Re + Lsum);
%             newQe = Qe(freeSigma,freeSigma);
%             newRL = RL(freeSigma) - newQe(freeSigma,lockedSigma)*RL(lockedSigma);
%
%             beta = zeros(obj.nStressDofs,1);



%             allSigma = 1:obj.nStressDofs;
%             freeSigma   = setdiff(allSigma, lockedSigma)
%
%             sigmaBc = [lockedSigma, lockedSigma*0];
%
%             ftemp = Re - LCAsum;
%             Qtemp = zeros(size(Qe));
%             keyboard;
%             Qtemp(freeSigma,freeSigma) = Qe(freeSigma,freeSigma)\(ftemp(freeSigma,:));
%
