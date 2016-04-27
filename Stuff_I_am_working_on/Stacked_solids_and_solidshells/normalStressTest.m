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
    
    eq = [0 0 0]';
    [Ke,fe] = soli8e(mesh.ex(:,elIndex)', mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', 2,elprop.D, eq);
%     [Ke,fe] = solid8(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex),elprop.D, eq,[3,3,3]);
    
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
end

%Add traction forces
f  = f + ftrac;

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

%Calulate stresses
ielitr = 1; figure; calfemzz = linspace(-1,1, mesh.nelz*2);
for iel=   6:11:(mesh.nelx*mesh.nelz)
    
    zz = [-1 0 1];
    for iz = 1:length(zz)
%     [strain(:,iz), stress(:,iz), xyzpoint(:,iz) ] = solid8Stress(mesh.ex(:,iel),mesh.ey(:,iel),mesh.ez(:,iel),ed(:,iel),elprop.D, [0,0,zz(iz)]);
%     [estmp ] = soli8s(mesh.ex(:,iel)',mesh.ey(:,iel)',mesh.ez(:,iel)',2,elprop.D,ed(:,iel)');
%     stress = [stress, estmp([1,5],:)'];
    end
    
    %---CALFEM
    xyzpoint = []; stress  =[];
    [estmp ] = soli8s(mesh.ex(:,iel)',mesh.ey(:,iel)',mesh.ez(:,iel)',2,elprop.D,ed(:,iel)');
    stress = [estmp([1,5],:)'];
    xyzpoint(3,:) = calfemzz( (1:2)  + (2*(ielitr-1)) )
    %----
    
    scmp = 3;
    plot(stress(scmp,:), xyzpoint(3,:)); hold on;
    
    ielitr = ielitr +1;
end




