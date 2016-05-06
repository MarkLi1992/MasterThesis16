function [ K,f, Ke, fe, He, Le ] = solid8eas(ex,ey,ez,D, M, eq,para)

%Integration rule
ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;

nEnhDofs = size(M(0,0,0),2);

%Init vectors
Ke = zeros(24,24);
He = zeros(nEnhDofs,nEnhDofs);
Le = zeros(nEnhDofs,24);
fe = zeros(24,1);

[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex', ey', ez');
T0 = transMat( JT0' );

%Start loop
for gp = ir.gps
   
    %Corods
    xi = gp.local_coords(1);eta = gp.local_coords(2);zeta = gp.local_coords(3);
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);
    N = solid8NMatrix(Nxieta);
    
    %Bmatrix
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex',ey',ez');
    B = solid8Bmatrix(dNdx);
    
    %Enanced part
    Mtemp = M(xi,eta,zeta);
    Mi = detJ0/detJ * T0*Mtemp;
    
    %Integrate
    Ke = Ke + B'*D*B * (detJ*gp.weight);
    He = He + Mi'*D*Mi * (detJ*gp.weight);
    Le = Le + Mi'*D*B * (detJ*gp.weight);

    fe = fe + N'*eq * (detJ*gp.weight);
    
end

K = Ke - Le'*inv(He)*Le;
f = fe;

end
