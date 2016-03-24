function [stresses, zcoords, xCoord] = solid8layered_postprocess(filename)

% load('SolidShell_x100_y50_z15.mat')
% load('Solid_x200_y100_z9.mat')
% load('SolidShell_x26_y13_z9.mat')
% filename = 'reference_solution_konsol/Solid_x200_y1_z10.mat';
load(filename)
 
%Element disp
ed = a(mesh.edof);

%Draw if you want to
iwanttotdraw = 0;
if(iwanttotdraw == 1)
    sfac = 100;
    exd = mesh.ex + ed(1:3:end,:)*sfac;
    eyd = mesh.ey + ed(2:3:end,:)*sfac;
    ezd = mesh.ez + ed(3:3:end,:)*sfac;
    
    fprintf('Drawing \n');
    figure(3);
    solid8draw(exd,eyd,ezd); hold on;
    view(3)
    axis equal
end

%Find element at coordinate
plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2])

%Get the stresses through the thickness at some points in a element
stressesAtPoints = [-1, 1 1 -1 0; -1 -1 1 1 0]*(1/sqrt(3));
[stresses, zcoords, xCoord] = el(plotEl).computeStressThroughThickness(a(mesh.edof(:,plotEl)), stressesAtPoints);
% [stresses, zcoords, xCoord] = el(1).computeStressThroughThickness(a(mesh.edof(:,plotEl)), stressesAtPoints);

%Plot all da shitz
stressComp = 5;
figure;
for ip=1:size(stressesAtPoints,2)
    subplot(2,3,ip)
    plot(stresses(ip).stress(stressComp,:), zcoords)
    title(sprintf('x=%.3f, y=%.3f', stressesAtPoints(1,ip), stressesAtPoints(2,ip) ))
end

% save('stressEas_0_90','stresses','zcoords');

