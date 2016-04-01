clear variables 

%Set up propblem
problem = 'KonsolMedUtbredd' %MedUtbredd
[mesh, elprop, M, bc, ftrac] = setup_problem(problem);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

for elIndex = 1:mesh.nel
    
    [Ke, fe] = solid8(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, [0,0,0]', [2,2,2]');
    
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
end

%Add traction forces
% f  = ftrac;

%Construct K-matrix
K = sparse(rows,cols,data);
clear rows cols data

%Solve equation of system
a = solveq(K,f,bc);
maxabs_a = max(abs(a));

ed = a(mesh.edof);
%Plot
exd = mesh.ex + ed(1:3:end,:);
eyd = mesh.ey + ed(2:3:end,:);
ezd = mesh.ez + ed(3:3:end,:);

figure;
solid8draw(exd,eyd,ezd); hold on;
view(3)
axis equal


%% LSF
sedof = zeros(8*6,mesh.nel);
Q = zeros(mesh.nno*6); L = zeros(mesh.nno*6, mesh.nno*3);
sdof = reshape(1:mesh.nno*6,6,mesh.nno)';
for i=1:mesh.nel
    m = sdof(mesh.nomesh(:,i),:)'; m = m(:);
    m2 = mesh.edof(:,i);
    sedof(:,i) = m;
    [ Qe, Le ] = LSFstress(mesh.ex(:,i), mesh.ey(:,i), mesh.ez(:,i), elprop.D);
    Q(m,m) = Qe; 
    L(m,m2) = Le;
end

as = Q\(L*a);

es = as(sedof);

figure
for i=1:mesh.nel
    
    cex = mesh.ex(:,i);
    ces = es(:,i);
%     ssx = ces([25,31]);
    ssx = ces([1,7]);
    xxx = cex([1,2]);
    
    plot(xxx,ssx); hold on;
    
end










