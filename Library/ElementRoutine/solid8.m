function [ Ke,fe ] = solid8(ex,ey,ez,D,eq,para)

ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;

%Init vectors
nDispDofs = 8*3;
Ke = zeros(nDispDofs,nDispDofs);
fe     = zeros(nDispDofs,1);

%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);
    N = solid8NMatrix(Nxieta);
    
    %Bmatrix
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex',ey',ez');
    B = solid8Bmatrix(dNdx);
    
    dV = detJ*gp.weight;
    Ke = Ke + B'*D*B * dV;
    fe = fe + N'*eq *  dV;
    
end

end