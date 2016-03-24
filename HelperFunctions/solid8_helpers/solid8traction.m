function [ f ] = solid8traction(ex,ey,t,ngp_)
%Ex and Ey are nodal coordinates for element
%t is the traction in form om a Anonomous function @(x,y)
%ir is the integration rule

if ~exist('ngp_')
    ngp_ = 2;
end

interp = InterpolatorX2Y2;
ir     = IntegrationRule2D();
ir.setupCubeRule(ngp_,ngp_);

f = zeros(12,1);

for gp = ir.gps
    
    lcoord = gp.local_coords;

    %Shape functions
    [Nxieta]=interp.eval_N(lcoord);
    
    N = [Nxieta(1)*eye(3), Nxieta(2)*eye(3),...
        Nxieta(3)*eye(3), Nxieta(4)*eye(3)];
    
    %Jacobian transpose (but ignore middle point)
    [~, detJ] = interp.eval_dNdx(lcoord, ex,ey);
    
    %Global x and y coord
    gx = Nxieta*ex';
    gy = Nxieta*ey';
    
    %Integrate
    dA = gp.weight * detJ;
    f = f + N'*t(gx,gy) * dA;
end


end


