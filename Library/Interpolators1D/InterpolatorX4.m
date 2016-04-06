classdef InterpolatorX4 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            z = coord(1);
            N = [ -(z/2 - 1/2)*((3*z)/2 + 1/2)*((3*z)/4 - 1/4), ((3*z)/2 - 1/2)*((3*z)/2 + 3/2)*((3*z)/4 - 3/4), -((3*z)/2 + 1/2)*((3*z)/2 - 3/2)*((3*z)/4 + 3/4), (z/2 + 1/2)*((3*z)/2 - 1/2)*((3*z)/4 + 1/4)];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, coord, ez)
            z = coord(1);
            dNdxi = [ - (3*(z/2 - 1/2)*((3*z)/2 + 1/2))/4 - (3*(z/2 - 1/2)*((3*z)/4 - 1/4))/2 - (((3*z)/2 + 1/2)*((3*z)/4 - 1/4))/2, (3*((3*z)/2 - 1/2)*((3*z)/2 + 3/2))/4 + (3*((3*z)/2 - 1/2)*((3*z)/4 - 3/4))/2 + (3*((3*z)/2 + 3/2)*((3*z)/4 - 3/4))/2, - (3*((3*z)/2 + 1/2)*((3*z)/2 - 3/2))/4 - (3*((3*z)/2 + 1/2)*((3*z)/4 + 3/4))/2 - (3*((3*z)/2 - 3/2)*((3*z)/4 + 3/4))/2, (3*(z/2 + 1/2)*((3*z)/2 - 1/2))/4 + (3*(z/2 + 1/2)*((3*z)/4 + 1/4))/2 + (((3*z)/2 - 1/2)*((3*z)/4 + 1/4))/2];
            
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

