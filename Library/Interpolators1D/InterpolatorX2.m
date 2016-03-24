classdef InterpolatorX2 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            xi = coord(1);
            N=[(-xi+1)/2 ,  (xi+1)/2];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, coord, ex)
            xi = coord(1);
            
            dNdxi=[-0.5, 0.5];
            
            JT=dNdxi*ex';
            detJ=det(JT);
            if (detJ<0)
                error('Jacobian not invertable')
            end
            
            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
        
    end
    
end

