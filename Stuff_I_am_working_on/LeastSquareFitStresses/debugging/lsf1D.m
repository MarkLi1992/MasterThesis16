function [ Me,Le ] = lsf1D( ex, g )

interp = InterpolatorX2;
ir = IntegrationRule;
ir.setupLineRule(2);

Me = zeros(2);
Le = zeros(2,1);

for gp = ir.gps
   
    N = interp.eval_N(gp.local_coords);
    [~,detJ] = interp.eval_dNdx(gp.local_coords, ex);
    
    dx = detJ*gp.weight;
    Me = Me + N'*N    * dx;  
    Le = Le + N'*g(1) * dx;
    
end


end

