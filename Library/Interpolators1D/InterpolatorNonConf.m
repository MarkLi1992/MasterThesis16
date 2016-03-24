classdef InterpolatorNonConf < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            xi = coord(1); eta = coord(2); zeta = coord(3);
            N=[1-xi^2, 1-eta^2, 1 - zeta^2];
        end
        
        function [dNdx] = eval_dNdx(obj, coord, ex, ey, ez)
            xi = coord(1); eta = coord(2); zeta = coord(3);
            
            dNdxi=[-2*xi, 0 0;...
                   0     -2*eta 0;...
                   0  , 0,  -2*zeta];
            
            interp3D = InterpolatorX2Y2Z2;
            [~, ~, JT] = interp3D.eval_dNdx(coord, ex, ey, ez);
            
            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
        
    end
    
end

