
%Element disp
ed = a(mesh.edof);

%Eb
% Iy = mesh.ly*mesh.lz^3/12; eb_maxdisp = abs((-1000)*mesh.lx^4)/(8*100e9*Iy); 
% maxabs_a = max(abs(a));

%Draw if you want to
iwanttotdraw = 1;
if(iwanttotdraw == 1)
    sfac = 1;
    exd = mesh.ex + ed(1:3:end,:)*sfac;
    eyd = mesh.ey + ed(2:3:end,:)*sfac;
    ezd = mesh.ez + ed(3:3:end,:)*sfac;
    
    fprintf('Drawing \n');
    figure;
    solid8draw(exd,eyd,ezd); hold on;
    view(3)
    axis equal
end

%Find element at coordinate
% plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2])
plotEl = 50%coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/4, mesh.ly/4, mesh.lz/2])
% plotEl = 50%805%90%50;
%Get the stresses through the thickness at some points in a element
% stressesAtPoints = [-1, 1 1 -1 0; -1 -1 1 1 0]*(1/sqrt(3));
stressesAtPoints = [0;0]*(1/sqrt(3));
[stresses, zcoords, xCoord] = el(plotEl).computeStressThroughThickness(ed(:,plotEl), stressesAtPoints);


%Plot all da shitz
stresscmpname = {'xx','yy','zz','xy','xz','yz'};
figure;
sitr = 1;
for stressComp = [1 2 4 5 6 3]
for ip=1:size(stressesAtPoints,2)
    subplot(2,3,sitr); sitr = sitr+1;
    plot(stresses(ip).stress(stressComp,:), zcoords,'-')
%     title(sprintf('x=%.3f, y=%.3f', stressesAtPoints(1,ip), stressesAtPoints(2,ip) ))
    title(sprintf('%s - element %i', stresscmpname{stressComp}, plotEl))
end
end

%save('SSsA3E3_hs2_61x61x10_AR10_stresses','stresses', 'zcoords')

% %Von mises
% for iel = 1:mesh.nel
%     [stress, zz] = el(iel).computeStressAt(ed(:,iel),[0,0,0]');
%     vm(iel) = sum(stress.stress.*(diag([1 1 1 2 2 2])*stress.stress));
%     ey2d(:,iel) = mesh.ey([1 5 7 3],iel);
%     ez2d(:,iel) = mesh.ez([1 5 7 3],iel);
% end
% figure;
% fill(ey2d,ez2d,vm)

% save('s_e3_stresses','stresses','zcoords')

% save('stressEas_0_90','stresses','zcoords');

