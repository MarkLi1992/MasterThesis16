clear variables;
% close all;

%Mesh data
lx = 0.1; ly=0.01; lz = 0.002; %lx = 0.1; ly=0.01; lz = 0.02;
nelx = 50; nely=5; nelz=1;
ndofsno = 9; %Number of dofs per node
P = -100; %force 4000
[edof,coord,ex,ey,ez,dof,nel,ndofs,nno, sideElements] = cubeMesher2(lx,ly,lz,nelx,nely,nelz,ndofsno);
neldofs = 8*ndofsno;

%Material data
E = 210e9; nu = 0.3;
D = hooke(4,E,nu);

%Assemble
n = nel*(neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1; f=zeros(ndofs,1);
for elIndex = 1:nel
    % Compute element stiffness
    [Ke,fe] = canon8e(ex(:,elIndex)', ey(:,elIndex)', ez(:,elIndex)', D, [0,0,-1*0]', [2,2,2]);
    elDofs = edof(:,elIndex);
    % Assemble
    for j = 1:neldofs
        for k = 1:neldofs
            rows(nPassed) = elDofs(j);
            cols(nPassed) = elDofs(k);
            data(nPassed) = Ke(j,k);
            nPassed = nPassed + 1;
        end
    end
    f(elDofs) = f(elDofs) + fe;
end

%Create K
K = sparse(rows,cols,data);

%Boundary force at end
[bc] = cubeBC('Konsol', dof(:,7:9), sideElements);
f(dof(sideElements(5).nodes(:),9)) = P/length(dof(sideElements(5).nodes(:),9));

%Solve
a = solveq(K,f,bc);

%Element disp.
ed = a(edof);

%Eulerbern
maxDisp = max(abs(a(9:9:end)));
II = (ly*lz^3/12);
EBmaxdisp = P*lx^3/3/E/II; %last term is I

fprintf('EulerBernoulli: %.5f, SolidElement: %.5f \n',EBmaxdisp,maxDisp);

% % %Get elementDisp
auedof = edof;
for i = 1:6
auedof(1:(9-i+1):end,:) = [];
end
ed = a(auedof);

asedof = edof;
for i = 1:3
asedof(7:(9-i+1):end,:) = [];
end
es = a(asedof);

figure
stressDofs = dof(:,1:7)';stressDofs = stressDofs(:);
as = a(stressDofs);
for i=1:nel
   
    cex = ex(:,i);
    ces = es(:,i);
%     ssx = ces([25,31]);
    ssxr = ces([25,31]);
    ssxl = ces([37,43]);
%     ssx = (ssxr + ssxl)./2;
    xxx = cex([1,2]);
    
    plot(xxx,ssxr); hold on;
    plot(xxx,ssxl);
end


%Plot
exd = ex + ed(1:3:end,:);
eyd = ey + ed(2:3:end,:);
ezd = ez + ed(3:3:end,:);

figure;
solid8draw(exd,eyd,ezd); hold on;
view(3)
axis equal

