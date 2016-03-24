classdef Solid8NonConform < handle
    
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
        interpNonConf
        ir
        
        %Submatrices are stored in case of EAS is used, and alpha is needed
        submatrices;
        
        %EAS matrix
        M;
    end
    
    methods
        
        function obj = Solid8NonConform(ex,ey,ez, elprop)
            
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
            
            obj.interpNonConf = InterpolatorNonConf;
            
            obj.lamIr = LayeredIntegrationRule(elprop.nLam, 2,2,2);
            
            obj.elprop = elprop;
            
        end
       
        function [K,f] = computeKandf(obj)
            
            nEnhDofs = 3*3;

            %Init vectors
            Kuu = zeros(24,24);
            Kaa = zeros(nEnhDofs,nEnhDofs);
            Kau = zeros(nEnhDofs,24);
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
                    dNdx_nonconf = obj.interpNonConf.eval_dNdx(ecoords, obj.lx(:,ilay)', obj.ly(:,ilay)', obj.lz(:,ilay)');
                    G = solid8Bmatrix( dNdx_nonconf );
                    
                    %Integrate
                    Kuu = Kuu + B'*D*B * (detJ*gp.weight);
                    Kaa = Kaa + G'*D*G * (detJ*gp.weight);
                    Kau = Kau + G'*D*B * (detJ*gp.weight);
                    
                    eq = [0,0,0]';
                    fe = fe + N'*eq * (detJ*gp.weight);
                end
            end
            K = Kuu - Kau'*inv(Kaa)*Kau;
            f = fe;
            
            obj.submatrices.Kuu = Kuu;
            obj.submatrices.Kaa = Kaa;
            obj.submatrices.Kau = Kau;
            
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
        
        function [alpha] = computeAlpha(obj, u)
            %Clarification of the code below: -----> alpha = inv(He)*Le*a
            alpha = -inv(obj.submatrices.Kaa)*obj.submatrices.Kau*u;
        end
        
        function [stresses, outcoords] = computeStressAt(obj,a,coords)
            
            %Get current layer at coordinates in "coords"
            clay = obj.getLayerFromZCoord(coords(3),'local');

            %D-matrix for this layer
            cD = obj.elprop.Dmatrices(:,:,clay);
            
            %Global coordinates
            outcoords = obj.interp.eval_N(coords) * [obj.ex,obj.ey,obj.ez];
            cAlpha = obj.computeAlpha(a);
            stresses = solid8NonConfStress(obj.ex,obj.ey,obj.ez, a, cAlpha, cD, coords, obj.interpNonConf, obj.interp);
            
        end
        
    
        
        function [Tx, Ty, outcoords] = computeShearForce(obj,a,xycoord)
            
            %if xycoord not specified, set it the the middle
            if~(exist('xycoord','var')); xycoord = [0,0]; end;
            
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
                    lcoords = [xycoord(1) xycoord(2) gp.local_coords];
                    
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
        
        function [] = deflectionInElement(obj,a)
            
            %Loop over x-dir
            nx = 20;
            xx = linspace(-1,1,nx);
            alpha = obj.computeAlpha(a);
            
            for ix = 1:nx
                coords = [xx(ix) 0 0];
                [Nxieta] = obj.interp.eval_N(coords)
                N = solid8NMatrix( Nxieta );
                
                [Pxieta] = obj.interpNonConf.eval_N(coords)
                P = solid8NMatrix(Pxieta);
                
                disp(:,ix) = P*alpha + N*a;
            end
           
            figure
            plot(xx,disp(3,:));
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



function [stress] = solid8NonConfStress(ex,ey,ez,a,alpha,D,coord, interpNonConf, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Shape functions
[dNdx, ~] = interp.eval_dNdx(coord, ex', ey', ez');

%Bmatrix
B = solid8Bmatrix(dNdx);

%NonconformingPart
[dNdx_nc] = interpNonConf.eval_dNdx(coord,ex',ey',ez');
G = solid8Bmatrix(dNdx_nc);

stress = D*B*a + D*G*alpha;
end

