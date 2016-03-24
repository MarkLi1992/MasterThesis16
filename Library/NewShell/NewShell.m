classdef NewShell
    %NEWSHELL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ir 
    end
    
    methods
        function ojb = NewShell()
           obj.ir = IntegrationRule;
           obj.ir.setup_cube_rule(2,2,2);
        end
        
        function [K,f] = computeKandf()

            
        end
    end
    
end

%The standard solid-elemenet stiffness
function [ Ke,fe ] = elementStiffness(ex,ey,ez,D,eq, interp, ir)

%Init vectors
Ke = zeros(24,24);
fe = zeros(24,1);

%Start loop
for gp= ir.gps
    
    %Shape functions
    Nxieta = interp.eval_N(gp.local_coords);
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords, ex, ey, ez);
    [N,B] = solid8NandBmatrix(Nxieta,dNdx);
    
    %derivative
    Ke = Ke + B'*D*B * (detJ*gp.weight);
    fe = fe + N'*eq * (detJ*gp.weight);
    
end

end

