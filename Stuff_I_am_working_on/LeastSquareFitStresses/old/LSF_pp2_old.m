

%mesh and edof for stresses
[ips_sedof,~,~,~,~,~,~,~,~,ips_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,2, 3);
[tau_sedof,~,~,~,~,~,~,~,~,tau_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,3, 2);

ips_a = zeros(max(max(ips_sedof)),el(1).elprop.nLam);

%Get in-plane stresses for all layers and for all elements
for i=1:mesh.nel
    stresses = el(i).computeStressThroughLayers(ed(:,i), [0 0]);
    for ilay=1:el(1).elprop.nLam
        es_xx(:,i,ilay) = stresses(1,1:2,ilay)';
        es_yy(:,i,ilay) = stresses(2,1:2,ilay)';
        es_xy(:,i,ilay) = stresses(4,1:2,ilay)';
    end
end

%Smear
as_xx = smearToNodes_MultipleLayers(ips_mesh, es_xx, 2);
as_yy = smearToNodes_MultipleLayers(ips_mesh, es_yy, 2);
as_xy = smearToNodes_MultipleLayers(ips_mesh, es_xy, 2); 

%Combine
ips_a(1:3:end,:) = as_xx;ips_a(2:3:end,:) = as_yy;ips_a(3:3:end,:) = as_xy;

% %smear stresses to nodes
% as = smearInPlaneStressesToNodes(mesh, el, a);

%Order in-plane-stresses in special format
es_layers = zeros(24,mesh.nel,el(1).elprop.nLam);
for i=1:el(1).elprop.nLam
    temp = ips_a(:,i);
    es_layers(:,:,i) = temp(ips_sedof);
end

%Get postprocessed shear stresses
postEl = 5;
for iel=1:mesh.nel
    es_input = [];
    for ilay=1:el(1).elprop.nLam
        es_input = [es_input, es_layers(:,iel,ilay)];
    end
    [temp_xz(:,iel), temp_yz(:,iel)] = el(iel).ShearStressesPostProcess(es_input, 'shear'); 
end
es_xz(:,:,1) = temp_xz(1:3,:);% es_xz(:,:,2) = temp_xz(3:5,:);
es_yz(:,:,1) = temp_yz(1:3,:);%es_yz(:,:,2) = temp_yz(3:5,:);

as_xz = smearToNodes_MultipleLayers(tau_mesh, es_xz, 3)
as_yz = smearToNodes_MultipleLayers(tau_mesh, es_yz, 3)

tau_a = zeros(max(max(tau_sedof)), el(1).elprop.nLam);
tau_a(1:2:end,:) = as_xz; tau_a(2:2:end,:) = as_yz;

for iel=1:mesh.nel
    input = [];
    for ilay=1:el(1).elprop.nLam
        temp = tau_a(:,ilay);
        tau_es = temp(tau_sedof);
        input = [input, tau_es(:,iel)];
    end
    temp_zz(:,iel) = el(iel).NormalStressPostProcess(es_input); 
end




% return
% as = smearStrainsToNodes(mesh, el, a);
% es = as(sedof);
% el(postEl).postProcess2(es(:,postEl))



% sdim = 3;
% sedof = zeros(sdim*8,mesh.nel);
% sdof = reshape(1:mesh.nno*sdim,sdim,mesh.nno)';
% for i=1:mesh.nel
%     m = sdof(mesh.nomesh(:,i),:)'; m = m(:);
%     sedof(:,i) = m;
% end

