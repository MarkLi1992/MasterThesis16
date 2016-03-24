function [ Kout, fout] = solid27e(ex,ey,ez,D,eq,para)

ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interpQuad = InterpolatorX3Y3Z3;

%Init vectors
nDispDofs = 27*3;

Ke = zeros(nDispDofs,nDispDofs);
fe     = zeros(nDispDofs,1);

%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interpQuad.eval_N(gp.local_coords);
    N = solid8NMatrix(Nxieta);
    
    %Jacobian transpose
    [dNdx, detJ] = interpQuad.eval_dNdx(gp.local_coords,ex,ey,ez);
    
    %Bmatrix
    B = solid8Bmatrix(dNdx);
     
    %Integrate
    dV  = gp.weight * detJ;
    Ke = Ke + B'*D*B * dV;
    
    fe = fe + N'*eq * dV;
end


% Kout = [zeros(24,24),Ce'; Ce, -Se];
% fout = [fe; zeros(48,1)];
Kout = Ke;
fout = fe;




end
