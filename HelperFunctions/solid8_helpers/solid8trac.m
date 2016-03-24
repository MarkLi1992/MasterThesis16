function [ f ] = solid8trac(ex,ey,t,ir)

if ~exist('ir')
    ir = 2;
end

ngpy = ir;
ngpx = ir;

%Get gauess weights and points
[gpx,gwx]=gausQuadTable(ir);
[gpy,gwy]=gausQuadTable(ir);

f = zeros(12,1);

%Start loop
for ix=1:ngpx
    for iy=1:ngpy
        
        xi=gpx(ix);
        eta=gpy(iy);
        
        
        %Shape functions
        [Nxieta, Bxieta]=shapeFunk2dLin(xi,eta);

        N = [Nxieta(1)*eye(3), Nxieta(2)*eye(3),...
             Nxieta(3)*eye(3), Nxieta(4)*eye(3)];

        %Jacobian transpose (but ignore middle point)
        JT=Bxieta*[ex', ey'];
        detJ=det(JT);

        if (detJ<eps)
            error('Jacobian not invertable')
        end
        
        dA = gwx(ix) * gwy(iy) * detJ;
        f = f + N'*t * dA;
    end
end


end

function [N,B]=shapeFunk2dLin(xi,eta,nn)

nn = 4;
if nn==4
%Byter plats på 3 och 4
    N=1/4*[(xi-1)*(eta-1),... 
          -(xi+1)*(eta-1),...
       -(xi-1)*(eta+1),...
           (xi+1)*(eta+1)];

    B=1/4*[eta-1, -(eta-1),  -(eta+1), eta+1;
           xi-1, -(xi+1), -(xi-1), xi+1];

% 
%     N=1/4*[(xi-1)*(eta-1),... 
%           -(xi+1)*(eta-1),... 
%            (xi+1)*(eta+1),...
%           -(xi-1)*(eta+1)];
% 
%     B=1/4*[eta-1, -(eta-1), eta+1, -(eta+1);
%            xi-1, -(xi+1), xi+1, -(xi-1)];

else
    error('Only 4-noded quadrilateral implemented')
end

end
