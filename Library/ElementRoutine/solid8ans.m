function [ Ke,fe ] = solid8Ans(ex,ey,ez,D,eq,para)

%Integration rule
ir = IntegrationRule;
ir.setupCubeRule(para(1), para(2), para(3));

interp = InterpolatorX2Y2Z2;
interpNormal = InterpolatorX2Y2;

%Init vectors
Ke = zeros(24,24);
fe = zeros(24,1);

[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex', ey', ez');
T0 = transMat( JT0' );
% T0 = eye(6);

%Interpolation points for ANS
[B13a,B23b,B13c,B23d] = BAnsShear(interp,  ex', ey', ez', inv(T0));
[B33e,B33f,B33g,B33h] = BAnsNormal(interp, ex', ey' ,ez', inv(T0));

%Start loop
for gp = ir.gps
   
    %Corods
    xi = gp.local_coords(1);eta = gp.local_coords(2);zeta = gp.local_coords(3);
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);
    N = solid8NMatrix(Nxieta);
    
    %Bmatrix
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex',ey',ez');
    Btmp = inv(T0)*solid8Bmatrix(dNdx);

    %Interpolation matrix for ANS   
    B13tilde = 0.5*(1-eta)*B13a + 0.5*(1+eta)*B13c; 
    Btmp(5,:) = B13tilde;
    
    B23tilde = 0.5*(1-xi)*B23d + 0.5*(1+xi)*B23b; 
    Btmp(6,:) = B23tilde;
    
%     N33 = interpNormal.eval_N(gp.local_coords([1 2]));
%     B33tilde =  N33(1)*B33e + N33(2)*B33f + N33(3)*B33h + N33(4)*B33g;
%     Btmp(3,:) = B33tilde;

    B = T0*Btmp;
    
    dV = detJ*gp.weight;
    Ke = Ke + B'*D*B * dV;
    fe = fe + N'*eq *  dV;
    
end

end
