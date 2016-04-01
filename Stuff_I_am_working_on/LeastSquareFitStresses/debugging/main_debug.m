clear variables

load('onedim_debug_data.mat')
sedof = sedof(1:2,:)';
nel = 10; ndofs = nel+1;

ex = reshape(xvec,2,nel)';
g  = reshape(svec,2,nel)';

M = zeros(ndofs,ndofs); L = zeros(ndofs,1);
for i = 1:nel
   
    [ Me , Le ] = lsf1D( ex(i, 1:2), g(i,1) );
    
    m = sedof(i,:);
    M(m,m) = M(m,m) + Me;
    L(m) = L(m) + Le;
    
end

as = M\L;

figure
for i=1:nel
    m = sedof(i,:);
    xx = ex(i,:);
    ss = as(m);
    gg = g(i,:);
    
    plot(xx,gg);
    hold on;
    plot(xx,ss);
    
end
    
    
    
    