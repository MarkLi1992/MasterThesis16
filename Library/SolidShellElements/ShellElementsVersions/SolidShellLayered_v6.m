classdef SolidShellLayered_v6 < handle
    %SolidShellLayered Sold-shell element for a laminate with different ply
    %properties.
    
    %Version 6 of layered solid shell element.
    %One EAS-M-matrix for all layers, 
    %stress-mesh
    
    
    properties
        
        %Layered integration rule
        lir;% = IntegrationRule;
        
        %Interolationrules
        dispInterp;
        stressInterp;
        stressInterpOrder; %interpolation order  for each stresss
        
        submatrices;
        
        %Interpolation matrix for strains
        Mhat = -1;
        
        %Number of dofs
        nDispDofs = -1;
        nStrainDofs = -1;
        nStressDofs = -1;
        nLamStressDofs = -1
        
        %Element coords
        ez; ey; ex;
        
        %Number of laminas
        nLam;
        
        %The deterimatant of each layer (always concstant for each element)
        lay_detJ;
        
        %Mesh-parameters for the stress
        stressMesh;
        
        %Dont know if alpha and beta should be stored in the element, by I
        %do it here temporararly
        beta;
        alpha; %layer-wise
    end
    
    methods
        function obj = SolidShellLayered_v6(ngpx, ngpy, ngpz, ex, ey, ez, stressInterp, Mhat, elprop)
            
            %stressInterp should define the interpolation for the
            %transverse stresses
            if(length(stressInterp) ~=3); error('Wrong interpolations for the stresses'); end;
            
            %Number of layers
            obj.nLam = elprop.nLam;
            
            %Create integration rules
            obj.lir = LayeredIntegrationRule(elprop.nLam, ngpx, ngpy, ngpz);
            
            %Create Strain interpolation matrices
            obj.Mhat = Mhat;
            
            %Define interpolation for each stress component
            %NOTE: only for the transverse stress components
            obj.stressInterp = getInterpolator([2 2 2], [2 2 2],stressInterp );
            obj.stressInterpOrder = stressInterp;
            
            %Interpolation for displacement always the same
            obj.dispInterp = InterpolatorX2Y2Z2;

            %Create stress-mesh. Better would probobly to have this as an
            %input instead (for effiency)
            obj.stressMesh = stressMesher(obj.nLam, stressInterp, ex,ey,ez, elprop.int_coordsG);
            
            %Define number of dofs
            obj.nStrainDofs = size(Mhat(0,0,0),2);
            obj.nLamStressDofs = sum(stressInterp*2*2);
            obj.nStressDofs = max(max(obj.stressMesh.edof));
            obj.nDispDofs = 3*2*2*2; %*nnoz??
            
            %Save element coords
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
        end
        
        function [R, J] = computeRandJ(obj,a, eq, ftrac, elprop)
            %Init
            V  = 0;
            
            QQ = zeros(obj.nStressDofs,obj.nStressDofs);
            KK1 = zeros(obj.nDispDofs, obj.nDispDofs);
            KK2 = zeros(obj.nDispDofs, obj.nStressDofs);
            RR = zeros(obj.nStressDofs, obj.nDispDofs);
            LLi= zeros(obj.nStressDofs, obj.nStrainDofs);
            
            LCAsum = zeros(obj.nStressDofs, obj.nDispDofs);
            fe = zeros(obj.nDispDofs,1);
            
            Ai = zeros(obj.nStrainDofs, obj.nDispDofs);
            Ci = zeros(obj.nStrainDofs, obj.nStrainDofs);
            
            %--new
            [~, detJ0, JT0] = obj.dispInterp.eval_dNdx([0,0,0], obj.ex, obj.ey, obj.ez);
            T0 = transMat( JT0' );
            %--new
            
            %Calculate alpha values
            for ilay = 1:obj.nLam
                layD = elprop.Dmatrices(:,:,ilay);
                
                Li = zeros(obj.nLamStressDofs, obj.nStrainDofs);
                
                Qe = zeros(obj.nLamStressDofs, obj.nLamStressDofs);
                Ke1 = zeros(obj.nDispDofs, obj.nDispDofs);
                Ke2 = zeros(obj.nDispDofs, obj.nLamStressDofs);
                Re = zeros(obj.nLamStressDofs, obj.nDispDofs);
                
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
%                     keyboard
                    %Derivatives of n-vector
                    [dNdx, detJEl] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
                    
                    % Get detJ for the layer
                    layZ = [[1,1,1,1]*elprop.int_coordsG(ilay), [1,1,1,1]*elprop.int_coordsG(ilay+1)];
%                     layZ = [[1,1,1,1]*elprop.int_coordsG(ilay),...
%                             [1,1,1,1]*(elprop.int_coordsG(ilay) + elprop.int_coordsG(ilay+1))/2,...
%                             [1,1,1,1]*elprop.int_coordsG(ilay+1)];
                    [~, detJ] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, layZ);                  
                        
                    % Get N and B matrix used in FEM
                    [N, B] = solid8NandBmatrix(Nvec, dNdx);
                    Bx = B([1 2 4],:);
                    Bz = B([3,5,6],:);

                    % Enanced part, eval in layered-system
                    Mi = obj.Mhat(ecoords(1), ecoords(2), ecoords(3));
                    %--new
                    Mi = detJ0/detJEl * T0*Mi;
                    %--new
                    
                    DM = layD*Mi;
                    DMz = DM([3 5 6],:);
                    
                    %Parts of the Dmatrix
                    DB = layD*B;
                    DBx = DB([1 2 4],:); 
                    DBz = DB([3 5 6],:);
                    
                    %Stress part, eval in layered-system
                    P = obj.getStressPmatrix(lcoords);
                    
                    % Integration
                    dV = detJ * gp.weight;
                    
                    V = V + 1*dV;
                    Ai = Ai + Mi'*layD*B  * dV;
                    Ci = Ci + Mi'*layD*Mi  * dV;
                    Li = Li + P' *(DMz)  * dV;
                    Re = Re + P' *(DBz)  * dV;
                    Qe = Qe + P'*P * dV;
                    Ke1 = Ke1 + Bx'*DBx * dV;
                    Ke2 = Ke2 + Bz'*P     * dV;
                    
                    fe = fe + N'*eq   * dV;

                end

                m = obj.stressMesh.edof(:,ilay);
                QQ(m,m) = QQ(m,m) + Qe;
                RR(m,:) = RR(m,:) + Re;
                KK1 = KK1 + Ke1;
                KK2(:,m) = KK2(:,m) + Ke2;
                LLi(m,:) = LLi(m,:) + Li;
                
            end
            
            
            alpha = Ci\(-Ai*a);
            
            LCAsum = LCAsum + LLi*inv(Ci)*Ai;
            %add the traction vector to the total force vector.
            fext = fe + ftrac;
            
            %Want to solve: ...Qe\(Re + Lsum);..., but we have some bc:s
            %Bc on sigma
            [TxybDofs, TxytDofs]  = obj.getStressComponentTopAndBottonDofs(5);
            [TzybDofs, TzytDofs]  = obj.getStressComponentTopAndBottonDofs(6);
            [SzzbDofs, SzztDofs]  = obj.getStressComponentTopAndBottonDofs(3);
            lockedSigma = [TxybDofs, TxytDofs TzybDofs, TzytDofs, SzzbDofs, SzztDofs]';

%             lockedSigma = [obj.getStressComponentDofs(5,[1:4,  9:12 17:20]),...
%                            obj.getStressComponentDofs(3,[1:4,  9:12 17:20]),...
%                            obj.getStressComponentDofs(6,[1:4,  9:12 17:20])]';

            sigmaBc = [lockedSigma, lockedSigma*0];
            sigmaBc((end-3):end,2) = -1000;
            
            allSigma = 1:obj.nStressDofs;
            freeSigma = setdiff(allSigma, lockedSigma);
            freeSigma = allSigma;
            
            %Solve for beta and derivate of beta wrt a
            %Caluclate righthandside of equations
            dfbeta = RR - LCAsum;
            fbeta  = dfbeta*a;
            
            %Init beta and dbeta/da
            beta    = zeros(obj.nStressDofs,1);
            beta(sigmaBc(:,1)) = sigmaBc(:,2);
            dbetada = zeros(obj.nStressDofs,obj.nDispDofs);
            
            %Solve (dont forget the bc:s :-) )
%             beta   (freeSigma  ) = QQ(freeSigma,freeSigma)\(fbeta(freeSigma) - QQ(freeSigma,lockedSigma)*sigmaBc(:,2));
            beta   (freeSigma  ) = QQ(freeSigma,freeSigma)\(fbeta(freeSigma));
            dbetada(freeSigma,:) = QQ(freeSigma,freeSigma)\(dfbeta(freeSigma,:));
            
            R = KK1*a + KK2*beta;%; - fext;
            J = (KK1  + KK2*dbetada);
            
            %Save internal variables for element
            obj.beta = beta;
            obj.alpha = alpha;
            
        end
        
        function [stress, zcoords]     = computeStressThroughThickness(obj)
            
            stressItr = 1;
            for ilay = 1:obj.nLam
                zz = linspace(-1, 1, 10);
                for iz = 1:length(zz);
                    
                    lcoords = [-0,0,zz(iz)];
                    
                    %Calculate P-matrix
                    P = obj.getStressPmatrix(lcoords);
                    
                    %Get dofs for this layer
                    el_beta = obj.stressMesh.edof(:,ilay);
                    
                    stress(:,stressItr) = P*obj.beta(el_beta);
%                     zcoords(stressItr)  = obj.dispInterp.eval_N(lcoords)*obj.stressMesh.ez(:,ilay);
                    tempInterpolator = InterpolatorX2Y2Z2;
                    zcoords(stressItr)  = tempInterpolator.eval_N(lcoords)*obj.stressMesh.ez(:,ilay);
                    stressItr = stressItr + 1;
                end
                
            end
            
        end
        
        function [stress,zcoords]     = computeCompatibleStressThroughThickness(obj,a,elprop)
            
            stressItr = 1;
            for ilam = 1:obj.nLam
                D = elprop.Dmatrices(:,:,ilam);
                
                zz = linspace(-1,1,3);
                for iz = 1:length(zz);
                    %Gauss points
                    lcoord = [0,0, zz(iz)]';
                    ecoords = lcoord;
                    ecoords(3) = obj.lir.getElementGaussCoordinate(...
                        lcoord(3), elprop.int_coordsL(ilam:ilam+1));
                    
                    %Shape functions
                    [dNdx, ~] = obj.dispInterp.eval_dNdx(ecoords, obj.ex, obj.ey, obj.ez);
                    
                    %Bmatrix
                    B = solid8Bmatrix(dNdx);
                    
                    stress(:,stressItr) = D*B*a + D*obj.Mhat(ecoords(1),ecoords(2),ecoords(3))*obj.alpha;
%                     zcoords(stressItr)  = obj.dispInterp.eval_N(lcoord)*obj.stressMesh.ez(:,ilam);
                    tempInterpolator = InterpolatorX2Y2Z2;
                    zcoords(stressItr)  = tempInterpolator.eval_N(lcoord)*obj.stressMesh.ez(:,ilam);
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
            
            %Total number of dofs for each stress-component in the element
            remap_s = [-1 -1 1 -1 2 3];
            component = remap_s(component);

            ss = (obj.nLam + 1 + (obj.stressInterpOrder-2)*obj.nLam)*4;
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i);
                startFrom = startFrom + temp;
            end
            
            dofs = startFrom + nodes;
            
        end
       
        function [alldofs] = getStressComponentAll(obj, component)
            
            %Total number of dofs for each stress-component in the element
            ss = (obj.nLam + 1 + (obj.stressInterpOrder-2)*obj.nLam)*4;
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i);
                startFrom = startFrom + temp;
            end
            endsAt = startFrom  + ss(component);
            
            alldofs = (startFrom+1):endsAt;
            
            
        end
        
        function [bottomDofs, topDofs] = getStressComponentTopAndBottonDofs(obj, component)
            
            %"Remap" stress components since only transverse stresses are
            %interpolated
            remap_s = [-1 -1 1 -1 2 3];
            component = remap_s(component);

            %Get the index from where the stress component strats from
            ss = obj.stressInterpOrder*4;
            startFrom = 0;
            for i=1:(component-1)
                temp = ss(i);
                startFrom = startFrom + temp;
            end
            endsAt = startFrom  + ss(component);
            
            temp = obj.stressMesh.edof(:,1);
            bottomDofs = temp((startFrom+1):(startFrom+4))';
            
            temp = obj.stressMesh.edof(:,end);
            topDofs = temp((endsAt-3):endsAt)';
            
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
% ss = (obj.nLam + 1 + (obj.stressInterpOrder-2)*obj.nLam)*4;
%             startFrom = 0;
%             for i=1:(component-1)
%                 temp = ss(i);
%                 startFrom = startFrom + temp;
%             end
%             endsAt = startFrom  + ss(component);
%
