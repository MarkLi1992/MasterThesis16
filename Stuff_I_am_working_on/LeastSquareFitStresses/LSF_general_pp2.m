
%mesh and edof for stresses
[ips_sedof,~,~,~,~,~,~,~,~,ips_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,2, 6);
[tau_sedof,~,~,~,~,~,~,~,~,tau_mesh,~] = cubeMesher(1,1,1,mesh.nelx,mesh.nely,mesh.nelz, 2,2,3, 2);

%Number of laminates
nLam = el(1).elprop.nLam;

%smear stresses to nodes
ips_a = smearInPlaneStressesToNodes(mesh, el, a);

%surf plot the streeses
es_xx = ips_a(2:6:end,1);
es_xx = es_xx(mesh.nomesh);
es_plot    = es_xx(1:4);
figure
fill(mesh.ey([1 5 7 3],:),mesh.ez([1 5 7 3],:), es_xx([1 5 7 3],:));  axis equal
% fill(mesh.ex([1 2 4 3],:),mesh.ey([1 2 4 3],:), es_xx([1 2 4 3],:));

%Get postprocessed shear stresses
clear tau_xz tau_yz
es_xz = zeros(3,mesh.nel,nLam); es_yz = zeros(3,mesh.nel,nLam);
for iel=1:mesh.nel
    %Input the shear stress post processor, Rows: The stresses in the nodes, Cols: layers 
    input = zeros(size(ips_sedof,1), nLam);
    for ilay=1:nLam
        temp = ips_a(:,ilay);
        ips_es = temp(ips_sedof);
        input(:,ilay) =ips_es(:,iel);
    end
%     figure('name',sprintf('e%i',iel))
    [tau_xz(:,iel), tau_yz(:,iel)] = el(iel).ShearStressesPostProcess2(input);%, 'shear'); 
    for ilay=1:nLam
        ind = (1:3) + 2*(ilay-1);
        es_xz(:,iel,ilay) = tau_xz(ind,iel);
        es_yz(:,iel,ilay) = tau_yz(ind,iel);
    end
end

as_xz = smearToNodes_MultipleLayers(tau_mesh, es_xz, 3);
as_yz = smearToNodes_MultipleLayers(tau_mesh, es_yz, 3);

tau_a = zeros(max(max(tau_sedof)), el(1).elprop.nLam);
tau_a(1:2:end,:) = as_xz; tau_a(2:2:end,:) = as_yz;

%Post process normal stress in transverse direction
clear sig_zz
es_zz = zeros(4,mesh.nel,nLam);
for iel=1:mesh.nel
    input = [];
    for ilay=1:el(1).elprop.nLam
        temp = tau_a(:,ilay);
        tau_es = temp(tau_sedof);
        input = [input, tau_es(:,iel)];
    end
    sig_zz(:,iel) = el(iel).NormalStressPostProcess(input); 
    for ilay=1:nLam
        ind = (1:4) + 3*(ilay-1);
        es_zz(:,iel,ilay) = sig_zz(ind,iel);
    end
end


%%Plot
plotEl = 30%805%90%50% 203%[ceil((mesh.nelx*mesh.nely)/2)]%[6];
% plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/4, mesh.ly/4, mesh.lz/2])
% plotEl = plotEl;
figure
for iel = plotEl
    zpoints_local = linspace(-1,1,10);
    zpoints_global = linspace(0,mesh.lz,10);
    thicknessInterpolator = InterpolatorX2; tauInterp = InterpolatorX3; sigzzInterp = InterpolatorX4;
    zitr = 1;
    for ilay = 1:nLam
        for iz = 1:length(zpoints_local)
            plot_tauxz(zitr) = tauInterp.eval_N(zpoints_local(iz))*es_xz(:,iel,ilay);
            plot_tauyz(zitr) = tauInterp.eval_N(zpoints_local(iz))*es_yz(:,iel,ilay);
            plot_sigzz(zitr) = sigzzInterp.eval_N(zpoints_local(iz))*es_zz(:,iel,ilay);
            plot_zz(zitr) = thicknessInterpolator.eval_N(zpoints_local(iz))*el(iel).elprop.lamZCoordsG(:,ilay);
            zitr = zitr + 1;
            
        end
    end
    
    subplot(2,3,4)
    plot(plot_tauxz, plot_zz,'k-*');
    title(sprintf('sigmaxz - Element %i', iel))
    ylabel('Thickness'); xlabel('sigma')
    
    subplot(2,3,5)
    plot(plot_tauyz, plot_zz,'k-*');
    title(sprintf('sigmayz - Element %i',iel))
    ylabel('Thickness'); xlabel('sigma')
    
    subplot(2,3,6)
    plot(plot_sigzz, plot_zz,'k-*');
    title(sprintf('sigmazz - Element %i',iel))
    ylabel('Thickness'); xlabel('sigma')
end
%inplanestresses
% [stresses, zco] = el(plotEl).computeStressThroughThickness(ed(:,plotEl),[0; 0]);
[~, ~, stresses, zco] = el(plotEl).computeStressThroughLayers(ed(:,plotEl),[0 0]); zco = [zeros(2,size(zco,2)); zco]; stresses.stress = stresses;

subplot(2,3,3)
plot_sigxy = stresses.stress(4,:);
plot(plot_sigxy,zco(3,:),'k-*');
title(sprintf('sigmaxy - Element %i',iel))
ylabel('Thickness'); xlabel('sigma')

subplot(2,3,2)
plot_sigyy = stresses.stress(2,:);
plot(plot_sigyy,zco(3,:),'k-*');
title(sprintf('sigmayy - Element %i',iel))
ylabel('Thickness'); xlabel('sigma')

subplot(2,3,1)
plot_sigxx = stresses.stress(1,:);
plot(plot_sigxx,zco(3,:),'k-*');
title(sprintf('sigmaxx - Element %i',iel))
ylabel('Thickness'); xlabel('sigma')


% save('sr_e12_stresses','plot_tauxz','plot_sigzz','plot_zz','plot_sigxx', 'zco');




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

