clear variables;

% addpath(genpath('Library'))
% addpath(genpath('HelperFunctions'))

% problem = 'Konsol';
problem = 'KonsolMedUtbredd';
% problem = 'InspandPlatta';
% problem = 'HybridStress2';
% problem = 'MembraneUnsymmetric';

[mesh, elprop, M, bc, ftrac] = setup_problem(problem);
eq = [0 0 0]'; 

% Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

err = 1;
iteration = 1;
tol =1e-4;
a = zeros(mesh.ndofs,1); % solution vector
% load a;

while err > tol

%     f = zeros(mesh.ndofs,1); % load vector
    R = zeros(mesh.ndofs,1); % Residual vector
    nPassed = 1;
    
    wh = waitbar(0,'Assembling...');
    for elIndex = 1:mesh.nel
       
%          el(elIndex) = SolidShellLayered_v1(3,3,3, mesh.ex(:,elIndex)', ...
%             mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', [2 2 3,2,3,3], M, elprop);
%         
%         el(elIndex) = SolidShellLayered_v2(3,3,10, mesh.ex(:,elIndex)', ...
%             mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', [2 2 3,2,7,7], M, elprop);
          
%          el(elIndex) = SolidShellLayered_v3(3,3,3, mesh.ex(:,elIndex)', ...
%             mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', [2 2 3,2,3,3], M, elprop);


         el(elIndex) = SolidShellLayered_v6(3,3,10, mesh.ex(:,elIndex)', ...
            mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', [4,3,3], M, elprop);

        
        elDofs = mesh.edof(:,elIndex);
        ae = a(elDofs);
        fe_in = (ftrac(elDofs)); %zeros(size(ftrac(elDofs)));%ftrac(:,elIndex);
%         Re = el(elIndex).computeR(ae, eq, ftrac(elDofs), elprop);
%         Ke = el(elIndex).computeJ(ae,elprop);  
        [Re, Ke] = el(elIndex).computeRandJ(ae, eq, fe_in, elprop);

        % Assemble
        for j = 1:mesh.neldofs
            for k = 1:mesh.neldofs
                rows(nPassed) = elDofs(j);
                cols(nPassed) = elDofs(k);
                data(nPassed) = Ke(j,k);
                nPassed = nPassed + 1;
            end
        end
        %f(elDofs) = f(elDofs) + fe;
        R(elDofs) = R(elDofs) + Re;
        
        waitbar(elIndex/mesh.nel,wh,...
            sprintf('Assembling residual and tangent for element %d (iteration %i)',...
            elIndex, iteration))
    end
    delete(wh)
    
    
    %Boundary condition
    R = R - ftrac;
    R(bc(:,1)) = 0;
    err = norm(R);
    
    fprintf('Iteration %i: Error %.3e \n', iteration, err)
    %Create K
    K = sparse(rows,cols,data);

    %Solve for increment delta a
    da = solveq(K,-R,bc);
    a = a + da;
    
    iteration = iteration + 1;
end

%Element disp.
ed = a(mesh.edof);

%Plot
exd = mesh.ex + ed(1:3:end,:);
eyd = mesh.ey + ed(2:3:end,:);
ezd = mesh.ez + ed(3:3:end,:);

% figure(3);
% solid8draw(exd,eyd,ezd); hold on;
% view(2)
% axis equal

% close all;
PostProcess_nonlin




