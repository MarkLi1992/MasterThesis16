function maxabs_a = quarterCylinderMain(nnn)
% clear variables

problem = 'SR1';
% problem = 'S';

[mesh, elprop, M, bc, ftrac] = quartercylinder_setup(problem, nnn);

% figure(88);
% solid8draw(mesh.ex,mesh.ey,mesh.ez); hold on;
% tmpEl = mesh.sideElements(3).elements;
% solid8draw(mesh.ex(:,tmpEl),mesh.ey(:,tmpEl),mesh.ez(:,tmpEl),8,'r'); view([-90 0]);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

for elIndex = 1:mesh.nel
    
%     el(elIndex) = Solid8StressRecLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, M);
%     el(elIndex) = Solid8AnsEasSR(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, M);
    
    eq = [-1e8*0 0 0]';
%     [ Ke,fe ] = solid8ans(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, eq,[5,5,5]);
%     [ Ke,fe ] = solid8(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, eq,[5,5,5]);
    [ Ke,fe ] = solid8anseas(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D,M, eq,[5,5,5]);
%     [ Ke,fe ] = solid8eas(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D,M, eq,[5,5,5]);
%     [Ke, fe] = el(elIndex).computeKandf(eq);
    
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
    if( mod(elIndex, 1000) == 0)
     fprintf('Assembling for element %i',elIndex);
    end
end
f = f + ftrac;
K = sparse(rows,cols,data);

fprintf('Solving \n');
a = solveq(K,f,bc);
maxabs_a = max(abs(a));

%max defl at right boundary
tmpdof = unique(mesh.dof(mesh.sideElements(3).nodes(:),3));
max(abs(a(tmpdof)));
min(abs(a(tmpdof)));
ed = a(mesh.edof);

%Plot
% sfac = 1e0;
% exd = mesh.ex + ed(1:3:end,:)*sfac;
% eyd = mesh.ey + ed(2:3:end,:)*sfac;
% ezd = mesh.ez + ed(3:3:end,:)*sfac;
% 
% figure(88);
% solid8draw(exd,eyd,ezd); hold on;
% view([-90 0])
% axis equal


