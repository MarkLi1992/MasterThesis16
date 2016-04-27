function [ strain, stress, xyzpoint ] = solid8Stress(ex,ey,ez,a,D,point)

interp = InterpolatorX2Y2Z2;

%Derivatives of x and y
N = interp.eval_N(point);
% N = solid8NMatrix(Nxy);

xyzpoint(1) = N*ex; xyzpoint(2) = N*ey; xyzpoint(3) = N*ez;
Bxy = interp.eval_dNdx(point, ex',ey',ez');
B = solid8Bmatrix(Bxy);

%Static condenstation
strain = B*a;
stress = D*strain;
end
