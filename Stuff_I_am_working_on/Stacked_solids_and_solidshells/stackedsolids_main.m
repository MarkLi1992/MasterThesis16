clear variables 
%Set up propblem
problem = ['HybridStress2','_stacked']; %KonsolMedUtbredd   CurvedBeam   HybridStress2
fprintf('Meshing\n');
[mesh, elprop, M, bc, ftrac] = setup_problem(problem);
eq = [0 0 0]'; 

% figure; solid8draw(mesh.ex,mesh.ey,mesh.ez); %axis equal

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

%Get stiffnessmatrix onece
% el(1) = Solid8EasLayered(mesh.ex(:,1), mesh.ey(:,1), mesh.ez(:,1), elprop, mesh.nlamel, M);

fprintf('Assembling %i x %i elements\n',mesh.nel, mesh.nlamel);
for elIndex = 1:mesh.nel
    
    el(elIndex) = Solid8layered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, mesh.nlamel);
%     el(elIndex) = Solid8EasLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, mesh.nlamel, M);
%     el(elIndex) = Solid8AnsEasLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, mesh.nlamel, M);
    [Ke, fe] = el(elIndex).computeKandf(eq);
    
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
    f(elDofs) = f(elDofs) + fe;
%     f(elDofs) = f(elDofs) + fe + ftrac(:,elIndex);
%     fprintf('Assembling for element %i \n',elIndex);
end

%Add traction forces
f  = f + ftrac;

%Construct K-matrix
K = sparse(rows,cols,data);
clear rows cols data

%Solve equation of system
fprintf('Solving\n');
a = solveq(K,f,bc);

maxabs_a = max(abs(a))
fprintf('Completed\n');


