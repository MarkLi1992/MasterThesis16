classdef Solid8StressRecLayered < handle
    
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
        ir
        
        %Submatrices are stored in case of EAS is used, and alpha is needed
        submatrices;
        
        %EAS matrix
        M;
    end
    
    methods
        
        function obj = Solid8StressRecLayered(ex,ey,ez, elprop, M)
            
            for i=1:elprop.nLam
                obj.lx(:,i) = ex;
                obj.ly(:,i) = ey;
                obj.lz(:,i) = [[1 1 1 1]*elprop.int_coordsG(i) [1 1 1 1]*elprop.int_coordsG(i+1)]';
            end
            obj.ex = ex; obj.ey = ey; obj.ez = ez;
            
            %Always simplest ir and intetp
            obj.interp = InterpolatorX2Y2Z2;
            obj.ir = IntegrationRule;
            obj.ir.setupCubeRule(2,2,2);
            
            obj.lamIr = LayeredIntegrationRule(elprop.nLam, 2,2,2);
            
            obj.elprop = elprop;
            obj.M = M;
            
            nEnhDofs = size(M(0,0,0),2);
            obj.submatrices.Ke = zeros(24);
            obj.submatrices.He = zeros(nEnhDofs);
            obj.submatrices.Le = zeros(nEnhDofs,24);
        end
       
        function [K,f] = computeKandf(obj)
            
            nEnhDofs = size(obj.M(0,0,0),2);

            %Init vectors
            Ke = zeros(24,24);
            He = zeros(nEnhDofs,nEnhDofs);
            Le = zeros(nEnhDofs,24);
            fe = zeros(24,1);
            
            for ilay = 1:obj.elprop.nLam
                D = obj.elprop.Dmatrices(:,:,ilay);
                for gp = obj.lamIr.irs(ilay).gps
                    
                    lcoords = gp.local_coords;
                    ecoords = lcoords;
                    ecoords(3) = obj.lamIr.getElementGaussCoordinate(lcoords(3), obj.elprop.int_coordsL(ilay:ilay+1));
                    
                    %Shape functions
                    Nxieta = obj.interp.eval_N(ecoords);
                    [dNdx, ~] = obj.interp.eval_dNdx(ecoords, obj.ex', obj.ey', obj.ez');
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
        
        function [] = postProcess(obj, as)
            
            Taubc = [0,0]';
            for ilay = 1:obj.elprop.nLam
                
                zminus = obj.elprop.int_coordsG(ilay);
                zplus  = obj.elprop.int_coordsG(ilay+1);
                
                zpoints = linspace(zminus,zplus,6);
                for iz = 1:length(zpoints);
                    sig_xz(1:2, tauItr) = sigxx_x*(zpoints(iz)-zminus) + sigxy_x*(zpoints(iz)-zminus) + Taubc;
                    
                    zplot(1, tauItr) = zpoints(iz);
                    tauItr = tauItr+1;
                end
                
                Taubc = sig_xz(:,end);
            end
            
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

