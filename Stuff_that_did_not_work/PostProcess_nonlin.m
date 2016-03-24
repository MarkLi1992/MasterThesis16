%Max displacement
maxDisp = max(abs(a));

%Area moment of inertia
Iy = (mesh.ly*mesh.lz^3/12); 

%Analytical values from Euler-bernoulli
eb_maxdisp = @(P,ly,lx,E,Iy) abs((P)*lx^4)/(8*E*Iy); %P = [N/m]

% Compare deflections between FEM and EB
% eb_maxdisp = eb_maxdisp(P,mesh.ly,mesh.lx, E,Iy);
fprintf('EulerBernoulli: %.10f, SolidElement: %.10f\n',eb_maxdisp(50/mesh.lx,mesh.ly,mesh.lx, 100e9 ,Iy),maxDisp);

elementsToPlot = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2]); % 
for iel = elementsToPlot

    [stressComp,zcoordsComp]  = el(iel).computeCompatibleStressThroughThickness(a(mesh.edof(:,iel)), elprop);
    plotCompStress = stressComp(5,:);
    
    [stress,zcoords] = el(iel).computeStressThroughThickness();
    plotStress = stress(2,:);

    figure
    plot(plotStress,zcoords);%, plotZZglobal); 
    hold on;
    plot(plotCompStress,zcoordsComp,'*'); 
    legend('Stress','Compatible stress')
    keyboard;
    
end