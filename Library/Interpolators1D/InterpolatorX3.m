classdef InterpolatorX3 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            z = coord(1);
            N = [ z*(z/2 - 1/2), -(z - 1)*(z + 1), z*(z/2 + 1/2)];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, coord, ez)
            z = coord(1);
            dNdxi = [ z - 1/2, -2*z, z + 1/2];
            
            JT=dNdxi*ez';
            detJ=det(JT);
            if (detJ<0)
                error('Jacobian not invertable')
            end

            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
     
    end
    
end

