classdef InterpolatorX2Y2 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, coord)
            xi = coord(1);
            eta = coord(2);
            N=1/4*[(xi-1)*(eta-1),...
                -(xi+1)*(eta-1),...
                -(xi-1)*(eta+1),...
                (xi+1)*(eta+1)];
            
            B=1/4*[eta-1, -(eta-1),  -(eta+1), eta+1;
                xi-1, -(xi+1), -(xi-1), xi+1];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, coord, ex,ey)
            xi = coord(1);
            eta = coord(2);
            
            dNdxi=1/4*[eta-1, -(eta-1),  -(eta+1), eta+1;
                xi-1, -(xi+1), -(xi-1), xi+1];
            
            JT=dNdxi*[ex',ey'];
            detJ=det(JT);
            if (detJ<0)
                error('Jacobian not invertable')
            end
            
            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
        
    end
    
end

