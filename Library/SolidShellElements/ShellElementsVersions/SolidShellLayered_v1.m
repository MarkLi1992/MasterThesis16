classdef SolidShellLayered_v1 < handle
    
    %Version 1 of layered solid shell element.
    %One EAS-M-matrix for each layers
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
        function obj = SolidShellLayered_v1(ngpx, ngpy, ngpz, ex, ey, ez, stressInterp, Mhat, elprop)
            
            %Number of layers
            obj.nLam = elprop.nLam;
            
            %Create integration rules
            obj.lir = LayeredIntegrationRule(elprop.nLam, ngpx, ngpy, ngpz);
            
            %Create Strain interpolation matrices
            obj.Mhat = Mhat;
            
            %Define interpolation for each stress component
            obj.stressInterp = getInterpolator( [2 2 2 2 2 2], [2 2 2 2 2 2],stressInterp );
            obj.stressInterpOrder = stressInterp;
            
            %Define number of dofs
            obj.nLamStrainDofs = size(Mhat(0,0,0),2);
            obj.nStrainDofs = obj.nLamStrainDofs * obj.nLam;
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
            Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
            Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
            Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
            Qe = zeros(obj.nStressDofs, obj.nStressDofs);
            Ke = zeros(obj.nDispDofs, obj.nStressDofs);
            Re = zeros(obj.nStressDofs, obj.nDispDofs);
            Lsum = zeros(obj.nStressDofs,1);
            LCAsum = zeros(obj.nStressDofs,obj.nDispDofs);
            fe = zeros(obj.nDispDofs,1);
            
            %Calculate alpha values
            for ilay = 1:obj.nLam
                layD = elprop.Dmatrices(:,:,ilay);
                
                Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
                Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
                Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
                
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
                    
                    % Enanced part, eval in layered-system
                    Mi = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));
                    
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
                    Ai = Ai + Mi'*layD*B  * dV;
                    Ci = Ci + Mi'*layD*Mi  * dV;
                    Li = Li + P' *layD*Mi  * dV;
                    Re = Re + P' *layD*B * dV;
                    Qe = Qe + P'*P * dV;
                    Ke = Ke + B'*P * dV;
                    
                    fe = fe + N'*eq   * dV;
                    
                    %Caluclate compatible stress
                    %                     gp.internalParameters = layD*(B*a + Mi*alpha
                end
                alpha(:,ilay) = Ci\(-Ai*a);
                %                 Lsum = Lsum + Li*alpha(:,ilay);
                LCAsum = LCAsum + Li*inv(Ci)*Ai;
            end

            %add the traction vector to the total force vector.
            fext = fe + ftrac;
            
            %Want to solve: ...Qe\(Re + Lsum);..., but we have some bc:s
            %Bc on sigma
%             [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
            lockedSigma = [obj.getStressComponentTopAndBottonDofs(5), ...
                            obj.getStressComponentTopAndBottonDofs(6)]';

            
            sigmaBc = [lockedSigma, lockedSigma*0];

            allSigma = 1:obj.nStressDofs;
            freeSigma = setdiff(allSigma, lockedSigma);
            
            %Solve for beta and derivate of beta wrt a
            %Caluclate righthandside of equations
            dfbeta = Re - LCAsum;
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
        
        function [stress, zcoords]     = computeStressThroughThickness(obj)
            
            stressItr = 1;
            zz = linspace(-1, 1, 10*obj.nLam);
            for iz = 1:length(zz);
                
                lcoords = [0,0,zz(iz)];
                
                %Calculate P-matrix
                P = obj.getStressPmatrix(lcoords);
                
                %Get dofs for this layer
                
                stress(:,stressItr) = P*obj.beta;
                zcoords(stressItr)  = obj.dispInterp.eval_N(lcoords)*obj.ez';
                stressItr = stressItr + 1;
            end
            
        end
        
        function [stress,zcoords]     = computeCompatibleStressThroughThickness(obj,a,elprop)
            
            stressItr = 1;
            for ilam = 1:obj.nLam
                D = elprop.Dmatrices(:,:,ilam);
                
                zz = linspace(-1,1,10);
                for iz = 1:length(zz);
                    %Gauss points
                    lcoord = [0,0, zz(iz)]';
                    
                    %Shape functions
                    [dNdx, ~] = obj.dispInterp.eval_dNdx(lcoord, obj.ex, obj.ey, obj.ez);
                    
                    %Bmatrix
                    B = solid8Bmatrix(dNdx);
                    
                    stress(:,stressItr) = D*B*a + D*obj.Mhat(lcoord(1),lcoord(2),lcoord(3))*obj.alpha(:,ilam);
                    zcoords(stressItr)  = obj.dispInterp.eval_N(lcoord)*...
                        [[1,1,1,1]*elprop.int_coordsG(ilam) ,[1,1,1,1]*elprop.int_coordsG(ilam+1) ]';
                    stressItr = stressItr+1;
                end
            end
            
        end
        
        function Knum = computeNumTangent(obj, R0, ae, eq, etrac, elprop, delta)%Re, ae, eq, etrac, elprop, delta
            % delta - size of numerical perturbation
            neldofs = length(R0);
            Knum = zeros(neldofs);
            
            for i=1:neldofs
                apert = ae;
                apert(i) = apert(i) + delta;
                Rpert = obj.computeR(apert, eq, etrac, elprop);
                Knum(:,i) = (Rpert-R0) / delta;
            end
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
        
        function [dofs] = getStressComponentTopAndBottonDofs(obj, component)
            
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
            dofs = [bottomDofs, topDofs];
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

%   function R = computeR(obj,a, eq, ftrac, elprop)
%             %Init
%             V  = 0;
%             Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
%             Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
%             Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
%             Qe = zeros(obj.nStressDofs, obj.nStressDofs);
%             Ke = zeros(obj.nDispDofs, obj.nStressDofs);
%             Re = zeros(obj.nStressDofs, obj.nDispDofs);
%             Lsum = zeros(obj.nStressDofs,1);
%
%             fe = zeros(obj.nDispDofs,1);
%             %Calculate alpha values
%             for ilay = 1:obj.nLam
%                 layD = elprop.Dmatrices(:,:,ilay);
%
%                 Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
%                 Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
%                 Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
%
%                 for gp = obj.lir.irs(ilay).gps %obj.lir.irs(ilay).gps
%
%                     %Gauss coordinates in local-layer-system and
%                     %local-element-system (only need to map the z-coord
%                     %since the x and y-coords should be the same
%                     lcoords = gp.local_coords;
%                     ecoords = lcoords;
%                     ecoords(3) = obj.lir.getElementGaussCoordinate(...
%                         lcoords(3), elprop.int_coordsL(ilay:ilay+1));
%
%                     %N-vector, eval in element-system
%                     Nvec = obj.dispInterp.eval_N(ecoords);
%
%                     %Derivatives of n-vector
%                     [dNdx, detJEl] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
%
%                     % Get detJ for the layer
%                     layZ = [[1,1,1,1]*elprop.int_coordsG(ilay), [1,1,1,1]*elprop.int_coordsG(ilay+1)];
%                     [~, detJ] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, layZ);
%
%                     % Get N and B matrix used in FEM
%                     [N, B] = solid8NandBmatrix(Nvec, dNdx);
%
%                     % Enanced part, eval in layered-system
%                     Mi = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));
%
%                     %Stress part, eval in element-system
%                     P = [];
%                     % Varför inte är inte alla 6 spänningskomponenter med? /JB
%                     for is = 1:length(obj.stressInterp)
%                         NvecSigma = obj.stressInterp{is}.eval_N(ecoords);
%                         % TODO: Mycket långsamt med blkdiag! /JB
%                         P = blkdiag(P,NvecSigma);
%                     end
%
%                     % Integration
%                     dV = detJ * gp.weight;
%
%                     V = V + 1*dV;
%                     Ai = Ai + Mi'*layD*B  * dV;
%                     Ci = Ci + Mi'*layD*Mi  * dV;
%                     Li = Li + P' *layD*Mi  * dV;
%                     Re = Re + P' *layD*B * dV;
%                     Qe = Qe + P'*P * dV;
%                     Ke = Ke + B'*P * dV;
%
%                     fe = fe + N'*eq   * dV;
%
%                     %Caluclate compatible stress
% %                     gp.internalParameters = layD*(B*a + Mi*alpha
%                 end
%                 alpha(:,ilay) = Ci\(-Ai*a);
%                 Lsum = Lsum + Li*alpha(:,ilay);
%             end
%             %             keyboard;
%
%             %add the traction vector to the total force vector.
%             fext = fe + ftrac ;
%
%             %Want to solve: ...Qe\(Re + Lsum);..., but we have some bc:s
%             %Bc on sigma
%             [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
%             lockedSigma = [TxybDofs, TxytDofs]';
%             sigmaBc = [lockedSigma, lockedSigma*0];
%
%             %solve for beta
%             %             beta = solveq(Qe, (Re + Lsum));
%             beta = solveq(Qe,(Re*a + Lsum),sigmaBc);
%             R = Ke*beta - fext;
%
%             %Save internal variables for element
%             obj.beta = beta;
%             obj.alpha = alpha;
%
%         end
%
%         function J = computeJ(obj,a,elprop)
%
%             %Init
%             V  = 0;
%             Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
%             Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
%             Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
%             Qe = zeros(obj.nStressDofs, obj.nStressDofs);
%             Ke = zeros(obj.nDispDofs, obj.nStressDofs);
%             Re = zeros(obj.nStressDofs, obj.nDispDofs);
%             LCAsum = zeros(obj.nStressDofs,obj.nDispDofs);
%
%             fe = zeros(obj.nDispDofs,1);
%             %Calculate alpha values
%             for ilay = 1:obj.nLam
%                 layD = elprop.Dmatrices(:,:,ilay);
%
%                 Ai = zeros(obj.nLamStrainDofs, obj.nDispDofs);
%                 Ci = zeros(obj.nLamStrainDofs, obj.nLamStrainDofs);
%                 Li = zeros(obj.nStressDofs, obj.nLamStrainDofs);
%
%                 for gp = obj.lir.irs(ilay).gps %obj.lir.irs(ilay).gps
%
%                     %Gauss coordinates in local-layer-system and
%                     %local-element-system (only need to map the z-coord
%                     %since the x and y-coords should be the same
%                     lcoords = gp.local_coords;
%                     ecoords = lcoords;
%                     ecoords(3) = obj.lir.getElementGaussCoordinate(...
%                         lcoords(3), elprop.int_coordsL(ilay:ilay+1));
%
%                     %N-vector, eval in element-system
%                     Nvec = obj.dispInterp.eval_N(ecoords);
%
%                     %Derivatives of n-vector
%                     [dNdx, detJEl] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
%
%                     % Get detJ for the layer
%                     layZ = [[1,1,1,1]*elprop.int_coordsG(ilay), [1,1,1,1]*elprop.int_coordsG(ilay+1)];
%                     [~, detJ] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, layZ);
%
%                     % Get N and B matrix used in FEM
%                     [N, B] = solid8NandBmatrix(Nvec, dNdx);
%
%                     % Enanced part, eval in layered-system
%                     Mi = obj.Mhat(lcoords(1), lcoords(2), lcoords(3));
%
%                     %Stress part, eval in element-system
%                     P = [];
%                     % Varför inte är inte alla 6 spänningskomponenter med? /JB
%                     for is = 1:length(obj.stressInterp)
%                         NvecSigma = obj.stressInterp{is}.eval_N(ecoords);
%                         % TODO: Mycket långsamt med blkdiag! /JB
%                         P = blkdiag(P,NvecSigma);
%                     end
%
%                     % Integration
%                     dV = detJ * gp.weight;
%
%                     V = V + 1*dV;
%                     Ai = Ai + Mi'*layD*B * dV; %*a
%                     Ci = Ci + Mi'*layD*Mi  * dV;
%                     Li = Li + P' *layD*Mi  * dV;
%                     Re = Re + P' *layD*B * dV;  %*a
%                     Qe = Qe + P'*P * dV;
%                     Ke = Ke + B'*P * dV;
%                 end
%                 LCAsum = LCAsum + Li*inv(Ci)*Ai;
%             end
%
%             %Want to solve: ...Qe\(Re + Lsum);..., but we have some bc:s
%             %Bc on sigma
%             [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
%             lockedSigma = [TxybDofs, TxytDofs]';
%             allSigma = 1:obj.nStressDofs;
%             freeSigma = setdiff(allSigma, lockedSigma);
%
%
%             %             dbda = Qe\(Re - LCAsum);
%             %             dbda(lockedSigma,:) = 0;
%
%             ftemp = Re - LCAsum;
%             dbda = zeros(obj.nStressDofs,obj.nDispDofs);
%             %              keyboard;
%             dbda(freeSigma,:) = Qe(freeSigma,freeSigma)\ftemp(freeSigma,:);
%
%
%             J = Ke*dbda; %inv(Qe)*(Re - LCAsum);
%
%         end
%

% 
% 
%         function compStress = computeCompatibleStressInLaminate(obj, a, lcoord, elprop)
%             
%             %Chack to see if lcoord is inside element bounds
%             if (~all(lcoord <= 1 & lcoord >= -1))
%                 error('Not inside elementbounds')
%             end
%             
%             layD = elprop.Dmatrices(:,:,ilay);
%             
%             %Gauss coordinates in local-layer-system and
%             %local-element-system (only need to map the z-coord
%             %since the x and y-coords should be the same
%             ecoords = lcoord;
%             ecoords(3) = obj.lir.getElementGaussCoordinate(ilay,...
%                 lcoord(3), elprop.int_coordsL(ilay:ilay+1));
%             
%             %N-vector, eval in element-system
%             Nvec = obj.dispInterp.eval_N(ecoords);
%             
%             %Derivatives of n-vector
%             [dNdx, ~] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
%             
%             % Get N and B matrix used in FEM
%             [~, B] = solid8NandBmatrix(Nvec, dNdx);
%             
%             % Enanced part, eval in layered-system
%             Mi = obj.Mhat(lcoord(1), lcoord(2), lcoord(3));
%             
%             %CALCULATE DAT SHIT:
%             compStress = layD*(B*a + Mi*obj.alpha(:,ilay));
%             
%         end
%         
%         function compStress = computeCompatibleStress(obj, a, elementCoord, elprop)
%             
%             %Chack to see if lcoord is inside element bounds
%             if (~all(elementCoord <= 1 & elementCoord >= -1))
%                 error('Not inside elementbounds')
%             end
%             
%             %What laminate are we in?:
%             ilay = sum(elementCoord(3) >= elprop.int_coordsL);
%             if(ilay > elprop.nLam); ilay = elprop.nLam; end;
%             
%             layD = elprop.Dmatrices(:,:,ilay);
%             
%             %Gauss coordinates in local-layer-system and
%             %local-element-system (only need to map the z-coord
%             %since the x and y-coords should be the same
%             lamCoord = elementCoord;
%             lamCoord(3) = obj.lir.getLayerGaussCoordinate(elementCoord(3), elprop.int_coordsL(ilay:ilay+1));
%             
%             %N-vector, eval in element-system
%             Nvec = obj.dispInterp.eval_N(elementCoord);
%             
%             %Derivatives of n-vector
%             [dNdx, ~] = obj.dispInterp.eval_dNdx(elementCoord, obj.ex, obj.ey, obj.ez);
%             
%             % Get N and B matrix used in FEM
%             [~, B] = solid8NandBmatrix(Nvec, dNdx);
%             
%             % Enanced part, eval in layered-system
%             Mi = obj.Mhat(lamCoord(1), lamCoord(2), lamCoord(3));
%             
%             %CALCULATE DAT SHIT:
%             compStress = layD*(B*a + Mi*obj.alpha(:,ilay));
%             
%         end