function [ Kout, fout, Se, Ce] = canon8e(ex,ey,ez,D,eq,para)

ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;

%Init vectors
nDispDofs = 8*3;
nStressDofs = 8*6;
Ce = zeros(nStressDofs,nDispDofs);
Se = zeros(nStressDofs,nStressDofs);
fe     = zeros(nDispDofs,1);

Scomp = inv(D);
%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);%
    Nu = solid8NMatrix(Nxieta);
    Nsigma =  solid8NMatrix(Nxieta, 6);
    
    %Jacobian transpose (but ignore middle point)
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex,ey,ez);

    B = solid8Bmatrix(dNdx);
    
    
    %Integrate
    dV  = gp.weight * detJ;
    Ce = Ce + Nsigma'*B * dV;
    Se = Se + Nsigma'*Scomp*Nsigma * dV; 
    
    fe = fe + Nu'*eq * dV;
end


% Kout = [zeros(24,24),Ce'; Ce, -Se];
% fout = [fe; zeros(48,1)];
Kout = [ Se, -Ce; -Ce', zeros(24,24)];
fout = [ zeros(48,1); -fe];

order = [];
for i=1:8
   eli = [(1:6)+((i-1)*6), (nStressDofs + (1:3)+((i-1)*3) )];
   order = [order, eli];
end
Kout = Kout(order,order);
fout = fout(order);

%Static condenstation
% Kout = Ce'*inv(Se)*Ce;
% fout = fe;



end
