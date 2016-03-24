classdef InterpolatorX1Y1Z2 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            z = coord(1);
            N = [ 1/2 - z/2, z/2 + 1/2];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, lcoords, ez)
            dNdxi = [ -1/2, 1/2];
            
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

