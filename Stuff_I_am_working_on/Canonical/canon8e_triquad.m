function [ Kout, fout, Se, Ce] = canon8e_triquad(ex,ey,ez,D,eq,para)

ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interpLin = InterpolatorX2Y2Z2;
interpQuad = InterpolatorX3Y3Z3;

%Init vectors
nDispDofs = 27*3;
nStressDofs = 27*6;
Ce = zeros(nStressDofs,nDispDofs);
Se = zeros(nStressDofs,nStressDofs);
fe     = zeros(nDispDofs,1);


invD = inv(D);
%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interpQuad.eval_N(gp.local_coords);%
    Nu = solid8NMatrix(Nxieta);
    Nsigma =  solid8NMatrix(Nxieta, 6);
    
    %Jacobian transpose (but ignore middle point)
    [dNdx, detJ] = interpQuad.eval_dNdx(gp.local_coords,ex,ey,ez);
    
    %Bmatrix
    B = solid8Bmatrix(dNdx);
     
    %Integrate
    dV  = gp.weight * detJ;
    Ce = Ce + Nsigma'*B * dV;
    Se = Se + Nsigma'*invD*Nsigma * dV; 
    
    fe = fe + Nu'*eq * dV;
end


% Kout = [zeros(24,24),Ce'; Ce, -Se];
% fout = [fe; zeros(48,1)];
Kout = [ Se, -Ce; -Ce', zeros(nDispDofs,nDispDofs)];
fout = [ zeros(nStressDofs,1); -fe];

order = [];
for i=1:27
   eli = [(1:6)+((i-1)*6), (nStressDofs + (1:3)+((i-1)*3) )];
   order = [order, eli];
end
Kout = Kout(order,order);
fout = fout(order);

%Static condenstation
% Kout = Ce'*inv(Se)*Ce;
% fout = fe;



end
