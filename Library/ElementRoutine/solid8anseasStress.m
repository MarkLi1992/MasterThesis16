function [stress] = solid8anseasStress(ex,ey,ez,a,alpha,D,Mi,coord, interp)

%Gauss points
xi=coord(1);
eta=coord(2);
zeta=coord(3);

%Trandformation matrix
[~, detJ0, JT0] = interp.eval_dNdx([0,0,0], ex', ey', ez');
T0 = transMat( JT0' );

%Interpolation points for ANS
[B13a,B23b,B13c,B23d] = BAnsShear(interp,  ex', ey', ez', inv(T0));
[B33e,B33f,B33g,B33h] = BAnsNormal(interp, ex', ey' ,ez', inv(T0));

%Shape functions
[dNdx, detJ] = interp.eval_dNdx(coord, ex', ey', ez');

%Bmatrix
% B = solid8Bmatrix(dNdx);

%BMatrix ---ANS---  
Btmp = inv(T0)*solid8Bmatrix(dNdx);
B13tilde = 0.5*(1-eta)*B13a + 0.5*(1+eta)*B13c; 
Btmp(5,:) = B13tilde;

B23tilde = 0.5*(1-xi)*B23d + 0.5*(1+xi)*B23b; 
Btmp(6,:) = B23tilde;


interpNormal = InterpolatorX2Y2;
N33 = interpNormal.eval_N([xi eta]);
B33tilde =  N33(1)*B33e + N33(2)*B33f + N33(3)*B33h + N33(4)*B33g;
Btmp(3,:) = B33tilde;
B = T0*Btmp;

%EAS part
Mtmp = Mi(xi,eta,zeta);
M = detJ0/detJ * T0 * Mtmp;

[Jmid] = interp.eval_ContraBaseVectors([0,0,0], ex', ey', ez');
T = transMat( Jmid );

%Stresses
stress = D*B*a + D*M*alpha;
stress = (T^-1)*stress;
end