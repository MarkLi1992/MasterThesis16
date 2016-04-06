
%Least square fit stresses
%number of stress dofs (6=all stress components, 3 = in-plane components)
sdim = 1;
%Edof for stresses
sedof = zeros(sdim*8,mesh.nel);
%Int matrices
M = zeros(mesh.nno*sdim); L = zeros(mesh.nno*sdim,1);% mesh.nno*3);
%dof for stresses
sdof = reshape(1:mesh.nno*sdim,sdim,mesh.nno)';
%Assemble
for i=1:mesh.nel
    m = sdof(mesh.nomesh(:,i),:)'; m = m(:);
    sedof(:,i) = m;
    [ Me, Le ] = el(i).computeLSFmatrices(ed(:,i));
%        Le2 = gradest(@(x) LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, x, ed(:,i)), ones(8,1)*10);
%        Me2 = hessian(@(x) LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, x, ed(:,i)), ones(8,1)*10);
    M(m,m) = Me; 
%     L(m,m2) = Le;
    L(m) = Le;
end

as = M\L;
% myas = [-1 -0.5 0,...
%         -1 -0.5 0,...
%         1 0.5 0,...
%         1 0.5 0]'*1e6;
%   myas = [-1 -1 -1 -1 1 1 1 1]'*1.4988e5;
  
es = as(sedof);
myes = myas(sedof);
% es = myas(sedof);

%Check potential
Pot1 = 0; Pot2 = 0;
for i=1:mesh.nel
    pot = LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, es(:,i), ed(:,i));
    Pot1 = Pot1 + pot;   
    
    pot = LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, myes(:,i), ed(:,i));
    Pot2 = Pot2 + pot;
    
    myLe=gradest(@(x) LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, x, ed(:,i)), 2*es(:,i));
    Le = gradest(@(x) LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, x, ed(:,i)), es(:,i));
    Me = hessian(@(x) LSFpontential(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D, x, ed(:,i)), es(:,i));
end


%look at stresses
figure
svec = []; xvec = [];
for i=1:mesh.nel
    
   ca = ed(:,i);
   
    %Kinematic stress
   stress1 = el(i).computeStressAt(ca, [-1,-1,1]);
   stress2 = el(i).computeStressAt(ca, [1,-1,1]);
   ss = [stress1(1), stress2(1)];
   svec = [svec, ss];
   %LSF stress
   ces = es(:,i);
   lsfss = ces([1,4]);
   lsfss = ces([13,16]);
   %xcoord;
   cex = mesh.ex(:,i);
   xx = cex([1,2]);
   xvec = [xvec, xx'];
   plot(xx,ss); hold on;
   plot(xx,lsfss,'k');
end







