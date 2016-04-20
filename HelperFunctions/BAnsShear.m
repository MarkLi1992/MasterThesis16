function [B13a,B23b,B13c,B23d] = BAnsShear(interp, ex, ey ,ez, invT0)
    
    %Corner coords
    A = [0,-1,0]; B = [1,0,0]; C = [0,1,0]; D = [-1,0,0];
    
    %-
    Bxy = interp.eval_dNdx(A,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B13a = Bmat(5,:);
    
    %-
    Bxy = interp.eval_dNdx(C,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B13c = Bmat(5,:);
    
    %-
    Bxy = interp.eval_dNdx(B,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B23b = Bmat(6,:);
    
    %-
    Bxy = interp.eval_dNdx(D,ex,ey,ez);
    Bmat = invT0*solid8Bmatrix(Bxy);
    B23d = Bmat(6,:);
    
end
