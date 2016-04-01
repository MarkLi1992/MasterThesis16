classdef Solid8RoffesQuadLayered < handle
    
    properties
        elprop;
        lamIr
        ex; ey; ez;
        
        %Internal mesh properties
        lx;ly;lz;
        
        %Z coordinate for all elements, for example: [-1 -0.9 ... 1]
        %Used when computing stresses
        elementZcoords
        
        %Interpolator and Integrationrule for each element
        interp
        interpQuad
        ir
        
        nDispDofs;
        
        %Submatrices are stored in case of EAS is used, and alpha is needed
        submatrices;
        
        %EAS matrix
        M;
    end
    
    methods
        
        function obj = Solid8RoffesQuadLayered(ex,ey,ez, elprop, M)
            
            for i=1:elprop.nLam
                obj.lx(:,i) = ex;
                obj.ly(:,i) = ey;
                obj.lz(:,i) = [[1 1 1 1 1 1 1 1 1]*elprop.int_coordsG(i) [1 1 1 1 1 1 1 1 1]*elprop.int_coordsG(i+1)]';
            end
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
            %Always simplest ir and intetp
            obj.interp = InterpolatorX2Y2Z2;
            ovj.interpQuad = InterpolatorX3Y3Z3;
            obj.ir = IntegrationRule;
            obj.ir.setupCubeRule(2,2,2);
            
            obj.lamIr = LayeredIntegrationRule(elprop.nLam, 2,2,2);
            
            obj.elprop = elprop;
            obj.M = M;
            
            obj.nDispDofs = 27*3;
            nEnhDofs = size(M(0,0,0),2);
            obj.submatrices.Ke = zeros(obj.nDispDofs);
            obj.submatrices.He = zeros(nEnhDofs);
            obj.submatrices.Le = zeros(nEnhDofs, obj.nDispDofs);
        end
       
        function [K,f] = computeKandf(obj)
            
            nEnhDofs = size(obj.M(0,0,0),2);

            %Init vectors
            Ke = zeros( obj.nDispDofs, obj.nDispDofs);
            He = zeros(nEnhDofs,nEnhDofs);
            Le = zeros(nEnhDofs, obj.nDispDofs);
            fe = zeros( obj.nDispDofs,1);
            
            for ilay = 1:obj.elprop.nLam
                D = obj.elprop.Dmatrices(:,:,ilay);
                for gp = obj.lamIr.irs(ilay).gps
                    
                    lcoords = gp.local_coords;
                    ecoords = lcoords;
                    ecoords(3) = obj.lamIr.getElementGaussCoordinate(lcoords(3), obj.elprop.int_coordsL(ilay:ilay+1));
                    
                    %Shape functions
                    Nxieta = obj.interpQuad.eval_N(ecoords);
                    [dNdx, ~] = obj.interpQuad.eval_dNdx(ecoords, obj.ex', obj.ey', obj.ez');
                    [N,B] = solid8NandBmatrix(Nxieta,dNdx);

                    [~, detJ] = obj.interp.eval_dNdx(ecoords, obj.lx(:,ilay)', obj.ly(:,ilay)', obj.lz(:,ilay)');
                    %Enanced part
                    Mi = obj.M(ecoords(1),ecoords(2),ecoords(3));
                    
                    %Integrate
                    Ke = Ke + B'*D*B * (detJ*gp.weight);
                    He = He + Mi'*D*Mi * (detJ*gp.weight);
                    Le = Le + Mi'*D*B * (detJ*gp.weight);
                    
                    eq = [0,0,0]';
                    fe = fe + N'*eq * (detJ*gp.weight);
                end
            end
            K = Ke - Le'*inv(He)*Le;
            f = fe;
            
            obj.submatrices.Ke = Ke;
            obj.submatrices.Le = Le;
            obj.submatrices.He = He;
            
        end
           
        function [stresses, outcoords] = computeStressThroughThickness(obj,a, local_points)
            %input:
            % local_points = [x1 x2 ...
            %                 y1 y2 ...];
            
            npoints = size(local_points,2);
            
            %Loop over the points 
            for ip=1:npoints
                iitr = 1;
%                 keyboard;
                %Loop over thickness points;
                zz = linspace(-1,1,100);
                for iz = 1:length(zz)
                    currentCoord =  [ local_points(:,ip); zz(iz)];
                    [stresses(ip).stress(:,iitr), outcoords(:,iitr)] = obj.computeStressAt(a,currentCoord);
                    iitr = iitr + 1;
                end
                
            end
            
            
        end
        
        function [alpha] = computeAlpha(obj, a)
            %Clarification of the code below: -----> alpha = inv(He)*Le*a
            alpha = -inv(obj.submatrices.He)*obj.submatrices.Le*a;
        end
        
        function [stresses, outcoords] = computeStressAt(obj,a,coords)
            
            %Get current layer at coordinates in "coords"
            clay = obj.getLayerFromZCoord(coords(3),'local');

            %D-matrix for this layer
            cD = obj.elprop.Dmatrices(:,:,clay);
            
            %Global coordinates
            outcoords = obj.interp.eval_N(coords) * [obj.ex,obj.ey,obj.ez];
            cAlpha = obj.computeAlpha(a);
            stresses = solid8EasLayeredStress(obj.ex,obj.ey,obj.ez, a, cAlpha, cD, obj.M, coords, obj.interp);
            
        end
        
        function [Tx, Ty, outcoords] = computeShearForce(obj,a, XY)
            if ~(exist('XY','var')); XY = [0,0]; end;
            %Interpolation and integrationRule for integration
            interp1D = InterpolatorX2;
            ir1d = IntegrationRule;
            ngp_xi = 3;
            ir1d.setupLineRule( ngp_xi);
            
            %Start gauss-loop
            Tx = 0; Ty = 0;
            for ilay = 1:obj.elprop.nLam
                for gp = ir1d.gps

                    %Layer z-coordinates
                    lz = obj.elprop.int_coordsG([ilay, ilay+1]);
                    
                    %Gauss points in this layer coordinate-system
                    lcoords = [XY(1) XY(2) gp.local_coords];
                    
                    %Coordinates in the element-coordinate system
                    ecoords = lcoords;
                    ecoords(3) = obj.lamIr.getElementGaussCoordinate(lcoords(3), obj.elprop.int_coordsL(ilay:ilay+1));
                    
                    %Get detJ
                    [~, detJ] = interp1D.eval_dNdx(gp.local_coords, lz);
                    
                    %Get stresses
                    [stresses, outcoords] = obj.computeStressAt(a,ecoords');

                    %Integrate
                    dz = gp.weight * detJ;
                    Tx = Tx + stresses(5) * (dz);
                    Ty = Ty + stresses(6) * (dz);
                end
            end

        end
        
        function [] = plzPostProcces(obj,a)
            
            %Stiffness matrices in 2D
            QLT = obj.elprop.D([1 2 4],[1 2 4]);
            QLTtilde = eye(2);
            
            %ABD matrices
            [A,B,D,Atilde] = ABDAtilde(QLT,QLTtilde, obj.elprop.angles, obj.elprop.int_coordsG);
            Dhat = [A,B;B,D];
            
            %Compute shearforces
            [Tx, Ty] = computeShearForce(obj,a);
            invDhat = inv(Dhat);
            
            %Compute strains
            nhatx = invDhat*[0 0 0 Tx 0 0]';
            nhaty = invDhat*[0 0 0 0 Ty 0]';
            
            epsbar_x = nhatx(1:3);
            kappa_x = nhatx(4:6);
            epsbar_y = nhaty(1:3);
            kappa_y = nhaty(4:6);
            
            %Integration rule and interpolation for all layers
            ir1d = IntegrationRule;
            ir1d.setupLineRule(3);
            interp1D = InterpolatorX2;
            
            %Derivatives of episolon
            %             eps_x = @(z) nhatx(1:3) + z*nhatx(4:6);
            %             eps_y = @(z) nhaty(1:3) + z*nhaty(4:6);

            %Boolean matrices
            Bx = [1 0 0; 0 0 1]; By = [0 0 1; 0 1 0];
            
            TauFunk = @(z, zminus, zmid, QLT) Bx*QLT*( z*(epsbar_x + kappa_x*(z*0.5-zmid))   -    zminus*(epsbar_x + kappa_x*(zminus*0.5-zmid)) ) +...
                                              By*QLT*( z*(epsbar_y + kappa_y*(z*0.5-zmid))   -    zminus*(epsbar_y + kappa_y*(zminus*0.5-zmid)) );
            
%             TauFunk = @(z, zminus, QLT) Bx*QLT*z*(epsbar_x + 0.5*z*kappa_x) + By*QLT*z*(epsbar_y + 0.5*z*kappa_y) +...
%                            Bx*QLT*zminus*(epsbar_x + 0.5*zminus*kappa_x) + By*QLT*zminus*(epsbar_y + 0.5*zminus*kappa_y);
            
            %integrate
            Taubc = [0,0]';
            tauItr = 1;
            zmid = 0%(obj.elprop.int_coordsG(1) + obj.elprop.int_coordsG(end))/2;
%             zmid2 = (obj.elprop.int_coordsG(1) + obj.elprop.int_coordsG(end))/2;
            for ilay = 1:obj.elprop.nLam
               QLT = obj.elprop.Dmatrices([1 2 4],[1 2 4],ilay);
                       
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                
                zpoints = linspace(zminus,zplus,6);
                for iz = 1:length(zpoints);
                    Tau(1:2, tauItr) = TauFunk(zpoints(iz), zminus, zmid, QLT) + Taubc;
                    zplot(1, tauItr) = zpoints(iz);
                    tauItr = tauItr+1;
                end
                
                Taubc = Tau(:,end);
                
            end
            figure
            plot(-Tau(1,:), zplot);
            keyboard;
        end
        
        function [el] = getLayerFromZCoord(obj,zc,opt)
            
            if(strcmp(opt,'global'))
              zcoords = obj.elprop.int_coordsG;
            elseif(strcmp(opt,'local'))
              zcoords = obj.elprop.int_coordsL;
            end
            
            for iz = 1:length(zcoords)-1
               
                if(zc >= zcoords(iz) && zc <= zcoords(iz+1))
                    el = iz;
                    return;
                end
                
            end
            
        end
        
    end
    
end



function [stress] = solid8EasLayeredStress(ex,ey,ez,a,alpha,D,Mi,coord, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Shape functions
[dNdx, ~] = interp.eval_dNdx(coord, ex', ey', ez');

%Bmatrix
B = solid8Bmatrix(dNdx);

%EAS part
M = Mi(xi,eta,zeta);

stress = D*B*a + D*M*alpha;
end

