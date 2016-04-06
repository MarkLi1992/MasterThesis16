function [ Qe, Le ] = LSFstress(ex,ey,ez,D)

ir = IntegrationRule;
ir.setupCubeRule(2,2,2);

interp = InterpolatorX2Y2Z2;

%Init vectors
nDispDofs = 24;
nStressDofs = 8*6;
Le = zeros(nStressDofs,nDispDofs);
Qe = zeros(nStressDofs,nStressDofs);

%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);%
    P =  solid8NMatrix(Nxieta, 6);
    
    %Jacobian transpose (but ignore middle point)
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex',ey',ez');
    B = solid8Bmatrix(dNdx);

    %Integrate
    dV  = gp.weight * detJ;
    Le = Le + P'*D*B * dV;
    Qe = Qe + P'*P * dV; 

end


end

