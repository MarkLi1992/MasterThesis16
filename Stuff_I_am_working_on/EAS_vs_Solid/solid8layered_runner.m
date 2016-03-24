function [maxabs_a] = solid8layered_runner(problem, ElementType, thickness, nelx, nely, nlamel)

% problem = 'Konsol';
% problem = 'InspandPlatta';
% problem = 'KonsolMedUtbredd';
% problem = 'TestUdiscont';

%Set up propblem
disp('Meshing');
[mesh, elprop, M, bc, ftrac, nlamel] = setup_problem_solidComposite(problem, thickness, nelx, nely, nlamel);
eq = [0 0 0]'; 

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1;
% f=zeros(mesh.ndofs,1);

% %Comute element stiffness once and for all
% switch ElementType
%     case 'SolidShell'
%        el(1) = Solid8EasLayered(mesh.ex(:,1), mesh.ey(:,1), mesh.ez(:,1), elprop, nlamel, M);
%     case 'Solid'
%        el(1) = Solid8layered(mesh.ex(:,1), mesh.ey(:,1), mesh.ez(:,1), elprop, nlamel);
%     otherwise
%         error('haha')
% end
% 
% [Ke, fe] = el(1).computeKandf();

fprintf('Assembling %i x %i elements\n',mesh.nel, nlamel);
for elIndex = 1:mesh.nel
   
    %Comute element stiffness once and for all
    switch ElementType
        case 'SolidShell'
           el(elIndex) = Solid8EasLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, nlamel, M);
        case 'Solid'
           el(elIndex) = Solid8layered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, nlamel);
        otherwise
            error('haha')
    end

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
%     f(elDofs) = f(elDofs) + fe;
%     fprintf('Assembling for element %i \n',elIndex);
end

%Add traction forces
f  = ftrac; %+ f;

%Construct K-matrix
K = sparse(rows,cols,data);
clear rows cols data

%Solve equation of system
fprintf('Solving\n');
a = solveq(K,f,bc);

maxabs_a = max(abs(a));
%Save and quit
% save(sprintf('reference_konsol_sxx/%s_x%i_y%i_z%i',ElementType,nelx,nely,nlamel));
fprintf('Completed\n');


