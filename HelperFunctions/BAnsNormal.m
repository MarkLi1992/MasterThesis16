function [B33e,B33f,B33g,B33h] = BAnsNormal(interp, ex, ey ,ez, invT0)
    
    %Corner coords
    E = [-1,-1,0]; F = [1,-1,0]; G = [1,1,0]; H = [-1,1,0];
    
    %-
    Bxy = interp.eval_dNdx(E,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B33e = Bmat(3,:);
    
    %-
    Bxy = interp.eval_dNdx(F,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B33f = Bmat(3,:);
    
    %-
    Bxy = interp.eval_dNdx(G,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B33g = Bmat(3,:);
    
    %-
    Bxy = interp.eval_dNdx(H,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B33h = Bmat(3,:);
    
end