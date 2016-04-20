
clear variables
lx = 2; inner_radius = 0.1; thickness = 0.005;
nelr = 1; nelphi = 30; nelx = 20;
[edof,coord,ex,ey,ez,dof,nel,ndofs,nno, nomesh, sideElements] = cylinderMesher( lx, inner_radius, thickness, nelr, nelphi, nelx, 2, 2, 2, 3);
mesh.edof = edof;mesh.nel = nel;mesh.dof = dof;mesh.coord = coord;mesh.ex = ex;mesh.ey = ey;mesh.ez = ez;mesh.ndofs = ndofs;mesh.nomesh = nomesh;mesh.sideElements = sideElements;mesh.neldofs = size(mesh.edof,1);

E = 100e9; nu = 0;
P = 9000;

figure; solid8draw(ex,ey,ez); axis equal; view(3)

elprop = ElementProperties();
elprop.setup(thickness, 0, hooke(4,E, nu));
nlamel = 1;
M = getInterPolMatrix(1);

%Assemble
n = mesh.nel*(mesh.neldofs)^2;
rows = zeros(n,1);cols = zeros(n,1);data = zeros(n,1);
nPassed = 1; f=zeros(mesh.ndofs,1);

for elIndex = 1:mesh.nel
    
    el(elIndex) = Solid8EasLayered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, nlamel, M);
%      el(elIndex) = Solid8layered(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop, nlamel);
    eq = [0 0 0]';
    [Ke, fe] = el(elIndex).computeKandf(eq);
%      [Ke, fe] = solid8(mesh.ex(:,elIndex), mesh.ey(:,elIndex), mesh.ez(:,elIndex), elprop.D, eq, [2,2,2]);
    
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
    f(elDofs) = f(elDofs) + fe; %ftrac(:,elIndex);
%     fprintf('Assembling for element %i \n',elIndex);
end

K = sparse(rows,cols,data);

fno = mesh.dof(mesh.sideElements(2).nodes(:),3)'; fno = unique(fno(:));
f(fno) = f(fno) + P/length(fno);


bc = mesh.dof(mesh.sideElements(1).nodes(:),:)'; bc = unique(bc(:));
bc = [bc,bc*0];
a = solveq(K,f,bc);

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


%Check
Iy = pi*(thickness/2 + inner_radius)^3 * thickness;
eb_maxdisp =abs(P*lx^3/3/E/Iy)
maxabs_a = max(abs(a))

