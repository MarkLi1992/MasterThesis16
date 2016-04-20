function [ K,f, Ke, fe, He, Le ] = solid8anseas(ex,ey,ez,D, M, eq,para)

%Integration rule
ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;
interpNormal = InterpolatorX2Y2;

nEnhDofs = size(M(0,0,0),2);

%Init vectors
Ke = zeros(24,24);
He = zeros(nEnhDofs,nEnhDofs);
Le = zeros(nEnhDofs,24);
fe = zeros(24,1);

%Interpolation points for ANS
[B13a,B23b,B13c,B23d] = BAnsShear(interp,  ex', ey', ez');
[B33e,B33f,B33g,B33h] = BAnsNormal(interp, ex', ey' ,ez');

[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex', ey', ez');
T0 = transMat( JT0 );
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

    %Interpolation matrix for ANS   
%     B13tilde = 0.5*(1-eta)*B13a + 0.5*(1+eta)*B13c; 
%     B(5,:) = B13tilde;
%     
%     B23tilde = 0.5*(1-xi)*B23d + 0.5*(1+xi)*B23b; 
%     B(6,:) = B23tilde;
    
    N33 = interpNormal.eval_N(gp.local_coords([1 2]));
    B33tilde =  N33(1)*T0*B33e + N33(2)*T0*B33f + N33(3)*T0*B33h + N33(4)*T0*B33g;
%     B33tilde = 0.25*((1-xi)*(1-eta)*B33e  +  (1+xi)*(1-eta)*B33f   +   (1+xi)*(1+eta)*B33g   +   (1-xi)*(1+eta)*B33h);
    B(3,:) = B33tilde;
    
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
