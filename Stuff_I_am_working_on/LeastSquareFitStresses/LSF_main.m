clear variables 
%Set up propblem

problem = 'KonsolMedUtbredd' %MedUtbredd
[mesh, elprop, M, bc, ftrac] = setup_problem(problem);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

for elIndex = 1:mesh.nel
    
    el(elIndex) = Solid8LSF(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop);
    [Ke, fe] = el(elIndex).computeKandf();
    
    % Assemble
    elDofs = mesh.edof(:,elIndex);
    for j = 1:mesh.neldofs
        for k = 1:mesh.neldofs
            rows(nPassed) = elDofs(j);
            cols(nPassed) = elDofs(k);
            data(nPassed) = Ke(j,k);
            nPassed = nPassed + 1;
        end
    end
    f(elDofs) = f(elDofs) + fe + ftrac(:,elIndex);
%     fprintf('Assembling for element %i \n',elIndex);
end

%Add traction forces
% f  = ftrac;

%Construct K-matrix
K = sparse(rows,cols,data);
clear rows cols data

%Solve equation of system
fprintf('Solving\n');
a = solveq(K,f,bc);

maxabs_a = max(abs(a));
fprintf('Completed\n');

ed = a(mesh.edof);

%Plot
exd = mesh.ex + ed(1:3:end,:);
eyd = mesh.ey + ed(2:3:end,:);
ezd = mesh.ez + ed(3:3:end,:);

figure;
solid8draw(exd,eyd,ezd); hold on;
view(3)
axis equal

%Least square fit stresses
%Edof for stresses
sedof = zeros(24,mesh.nel);
%Int matrices
M = zeros(mesh.nno*3); L = zeros(mesh.nno*3,1);% mesh.nno*3);
%dof for stresses
sdof = reshape(1:mesh.nno*3,3,mesh.nno)';
%Assemble
for i=1:mesh.nel
    m = sdof(mesh.nomesh(:,i),:)'; m = m(:);
    m2 = mesh.edof(:,i);
    sedof(:,i) = m;
    [ Me, Le ] = el(i).computeLSFmatrices(ed(:,i));
    M(m,m) = Me; 
%     L(m,m2) = Le;
    L(m) = Le;
end

as = M\L;
es = as(sedof);
%look at stresses
figure
for i=1:mesh.nel
    
   ca = ed(:,i);
   
    %Kinematic stress
   stress1 = el(i).computeStressAt(ca, [-1,0,-1]);
   stress2 = el(i).computeStressAt(ca, [1,0,-1]);
   ss = [stress1(1), stress2(1)];
   
   %LSF stress
   ces = es(:,i);
   lsfss = ces([1,4]);
   
   %xcoord;
   cex = mesh.ex(:,i);
   xx = cex([1,2]);
   
   plot(xx,ss); hold on;
   plot(xx,lsfss,'k');
end







