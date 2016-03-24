
%Element disp
ed = a(mesh.edof);

%Draw if you want to
iwanttotdraw = 0;
if(iwanttotdraw == 1)
    sfac = 1;
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
[stresses, outCoords] = el(1).computeStressThroughThickness(a(mesh.edof(:,plotEl)), stressesAtPoints);
zcoords = outCoords(3,:);

%Plot all da stresses
stressComp = 1;
figure;
for ip=1:size(stressesAtPoints,2)
    subplot(2,3,ip)
    plot(stresses(ip).stress(stressComp,:), zcoords)
    title(sprintf('x=%.3f, y=%.3f', stressesAtPoints(1,ip), stressesAtPoints(2,ip) ))
end

%plot shear forces
figure;
for i=1:mesh.nel
    [Tx1, Ty1, outcoords1]  = el(i).computeShearForce(a(mesh.edof(:,i)),[-1 0])
    x1 = outcoords1(1);
    [Tx2, Ty2, outcoords2]  = el(i).computeShearForce(a(mesh.edof(:,i)), [1 0])
    x2 = outcoords2(1);
    
    plot([x1,x2],[Tx1,Tx2]); hold on;
end

%calculate element displacement, they should be quatratic
el(plotEl).deflectionInElement(a(mesh.edof(:,plotEl)));


