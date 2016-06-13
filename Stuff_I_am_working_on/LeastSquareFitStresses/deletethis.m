
%mesh and edof for stresses
[ips_sedof,~,~,~,~,~,~,~,~,ips_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,2, 6);
[tau_sedof,~,~,~,~,~,~,~,~,tau_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,3, 2);

%Number of laminates
nLam = el(1).elprop.nLam;

%smear stresses to nodes
ips_a = smearInPlaneStressesToNodes(mesh, el, a);

figure;
for iel = 1:mesh.nel
    [s1] = el(iel).computeStressThroughLayers(ed(:,iel), [0,-1]);
    [s2] = el(iel).computeStressThroughLayers(ed(:,iel), [0,1]);
    
    plot([0 1] + (iel), [s1(2,1) s2(2,1)]); hold on;
end

%surf plot the streeses
es_xx = ips_a(2:6:end,1);
es_xx = es_xx(mesh.nomesh);
es_plot    = es_xx(1:4);

plot(es_xx([1],:)); 

figure
fill(mesh.ey([1 5 7 3],:),mesh.ez([1 5 7 3],:), es_xx([1 5 7 3],:));  axis equal



