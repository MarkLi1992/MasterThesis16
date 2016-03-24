clear variables

problem = 'Konsol';
% problem = 'HybridStress2';

%Set up propblem
[mesh, elprop, M, bc, ftrac] = setup_problem(problem);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);
f = zeros(mesh.ndofs,1);
nPassed = 1;

% fprintf('Assembling %i x %i elements\n',mesh.nel, nlamel);
for elIndex = 1:mesh.nel
    el(elIndex) = Solid8RoffesLayered(mesh.ex(:,elIndex),mesh.ey(:,elIndex),mesh.ez(:,elIndex), elprop, M);
%     el(elIndex) = Solid8NonConform(mesh.ex(:,elIndex),mesh.ey(:,elIndex),mesh.ez(:,elIndex), elprop);
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
end

%Construct K-matrix
K = sparse(rows,cols,data);
clear rows cols data

%Solve equation of system
fprintf('Solving\n');
a = solveq(K,f,bc);

maxabs_a = max(abs(a));

plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2])

return;
for i=plotEl%1:mesh.nel
    [stresses, coords] = el(i).computeStressThroughThickness(a(mesh.edof(:,i)),[0;0]);
    figure
    plot(stresses.stress(1,:), coords(3,:))
    el(i).plzPostProcces(a(mesh.edof(:,i))); %computeShearForce(a(mesh.edof(:,i)))%
%    stress2 = el(i).computeStressAt(a(mesh.edof(:,i)),[0 0 1])
end



