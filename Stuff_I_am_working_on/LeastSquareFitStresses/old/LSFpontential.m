function [ pot ] = LSFpontential(ex,ey,ez,D, alpha, a)

ir = IntegrationRule;
ir.setupCubeRule(2,2,2);

interp = InterpolatorX2Y2Z2;

%Init vectors
sdim = 1;
nDispDofs = 24;
nStressDofs = 8*sdim;
pot = 0;

%Start loop
for gp = ir.gps
    
    %Shape functions
    [Nxieta]= interp.eval_N(gp.local_coords);%
    P =  solid8NMatrix(Nxieta, sdim);
    
    %Jacobian transpose (but ignore middle point)
    [dNdx, detJ] = interp.eval_dNdx(gp.local_coords,ex',ey',ez');
    B = solid8Bmatrix(dNdx);
    stress = D*B*a;
    
    %Integrate
    dV  = gp.weight * detJ;
    pot = pot + 0.5*(P*alpha - stress(1))'   *   (P*alpha - stress(1)) * dV;

end


end

