function [ Kout, fout, Qe, Le, Ce] = lsf8e(ex,ey,ez,D,eq,para)

ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;

%Init vectors
sdim = 6;
nDispDofs = 8*3;
nStressDofs = 8*sdim;
Ce = zeros(nDispDofs   ,nStressDofs);
Qe = zeros(nStressDofs ,nStressDofs);
Le = zeros(nStressDofs ,nDispDofs);

fe     = zeros(nDispDofs,1);

%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);%
    Nu = solid8NMatrix(Nxieta, 3);
    Nsigma =  solid8NMatrix(Nxieta, sdim);
    
    %Jacobian transpose (but ignore middle point)
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex,ey,ez);

    B = solid8Bmatrix(dNdx);
    
    
    %Integrate
    dV  = gp.weight * detJ;
    Ce = Ce + B'*Nsigma * dV;
    Qe = Qe + Nsigma'*Nsigma * dV; 
    Le = Le + Nsigma'*D*B * dV;
    
    fe = fe + Nu'*eq * dV;
end


Kout = [ Qe, -Le; Ce, zeros(nDispDofs,nDispDofs)];
fout = [ zeros(nStressDofs,1); fe];

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
