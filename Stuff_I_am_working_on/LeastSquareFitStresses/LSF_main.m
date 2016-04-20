clear variables 
addpath('smearers')
%Set up propblem

problem = 'KonsolMedUtbredd'; %KonsolMedUtbredd InspandPlatta
[mesh, elprop, M, bc, ftrac] = setup_problem(problem);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

for elIndex = 1:mesh.nel
    
%     el(elIndex) = Solid8LSF(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop);Solid8StressRecLayered
%     el(elIndex) = Solid8StressRecLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, M);
    eq = [0 0 0]';
%     [Ke, fe] = el(elIndex).computeKandf(eq);
%     [ Ke,fe ] = solid8ans(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, eq,[3,3,3]);
    [ Ke,fe ] = solid8anseas(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, M, eq,[3,3,3]);
    
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
fprintf('Solving ');
a = solveq(K,f,bc);

maxabs_a = max(abs(a));
fprintf('Completed\n');

ed = a(mesh.edof);

%Plot
sfac = 1;
exd = mesh.ex + ed(1:3:end,:)*sfac;
eyd = mesh.ey + ed(2:3:end,:)*sfac;
ezd = mesh.ez + ed(3:3:end,:)*sfac;

figure;
solid8draw(exd,eyd,ezd); hold on;
view(3)
axis equal
