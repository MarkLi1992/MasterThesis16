classdef Solid8LSF < handle
    
    properties
        elprop;
        lamIr
        ex; ey; ez;
        
        %Internal mesh properties
        lx;ly;lz;
        
        %Interpolator and Integrationrule for each element
        interp
        ir
        
        %Submatrices are stored in case of EAS is used, and alpha is needed
        submatrices;
    end
    
    methods
        
        function obj = Solid8LSF(ex,ey,ez,elprop)
            
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
            
        end
        
        function [Ke,fe] = computeKandf(obj)
            
            Ke = zeros(24,24); fe = zeros(24,1);
            for gp = obj.ir.gps
                
                %Shape functions
                [Nxieta]= obj.interp.eval_N(gp.local_coords);
                N = solid8NMatrix(Nxieta);
                
                %Bmatrix
                [dNdx, detJ] = obj.interp.eval_dNdx(gp.local_coords,obj.ex',obj.ey',obj.ez');
                B = solid8Bmatrix(dNdx);
                
                dV = detJ*gp.weight;
                Ke = Ke + B'*obj.elprop.D*B * dV;
                eq = [0,0,0]';
                fe = fe + N'*eq *  dV;
                
            end
            
            obj.submatrices.Ke = Ke;
            
        end
        
        function [Me,Le] = computeLSFmatrices(obj, a)
            
            Me = zeros(24,24);  Le = zeros(24,1);%,24);
            %Start loop
            for gp = obj.ir.gps
                
                %Shape functions
                [Nxieta]= obj.interp.eval_N(gp.local_coords);%
                P =  solid8NMatrix(Nxieta, 3);
                
                %Jacobian transpose (but ignore middle point)
                [dNdx, detJ] = obj.interp.eval_dNdx(gp.local_coords,obj.ex',obj.ey',obj.ez');
                B = solid8Bmatrix(dNdx);
                stress = obj.elprop.D*B*a;
                
                %Integrate
                dV  = gp.weight * detJ;
                Le = Le + P'*stress([1 2 4]) * dV;
                Me = Me + P'*P * dV;
                
            end
            
            obj.submatrices.Me = Me;
            obj.submatrices.Le = Le;
            
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
        
        
        function [stresses, outcoords] = computeStressAt(obj,a,coords)
            
            %Get current layer at coordinates in "coords"
            clay = obj.getLayerFromZCoord(coords(3),'local');

            %D-matrix for this layer
            cD = obj.elprop.Dmatrices(:,:,clay);
            
            %Global coordinates
            outcoords = obj.interp.eval_N(coords) * [obj.ex,obj.ey,obj.ez];
            stresses = solid8Stress(obj.ex,obj.ey,obj.ez, a, cD, coords, obj.interp);
            
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



function [stress] = solid8Stress(ex,ey,ez,a,D,coord, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Shape functions
[dNdx, ~] = interp.eval_dNdx(coord, ex', ey', ez');

%Bmatrix
B = solid8Bmatrix(dNdx);

stress = D*B*a;
end

