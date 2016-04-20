
%Element disp
ed = a(mesh.edof);

%Eb
Iy = mesh.ly*mesh.lz^3/12; eb_maxdisp = abs((-50)*mesh.lx^4)/(8*100e9*Iy); 
maxabs_a = max(abs(a));

%Draw if you want to
iwanttotdraw = 0;
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
plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2])
% plotEl = 203%805%90%50;
%Get the stresses through the thickness at some points in a element
% stressesAtPoints = [-1, 1 1 -1 0; -1 -1 1 1 0]*(1/sqrt(3));
stressesAtPoints = [0; 0]*(1/sqrt(3));
[stresses, zcoords, xCoord] = el(1).computeStressThroughThickness(ed(:,plotEl), stressesAtPoints);


%Plot all da shitz
stresscmpname = {'xx','yy','zz','xy','xz','yz'};
figure;
sitr = 1;
for stressComp = [1 2 4 5 6 3]
for ip=1:size(stressesAtPoints,2)
    subplot(2,3,sitr); sitr = sitr+1;
    plot(stresses(ip).stress(stressComp,:), zcoords)
%     title(sprintf('x=%.3f, y=%.3f', stressesAtPoints(1,ip), stressesAtPoints(2,ip) ))
    title(sprintf('sigma_{%s} - element %i', stresscmpname{stressComp}, plotEl))
end
end
% save('s_e3_stresses','stresses','zcoords')

% save('stressEas_0_90','stresses','zcoords');

