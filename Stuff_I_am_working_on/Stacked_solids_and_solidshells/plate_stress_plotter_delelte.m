
ed = a(mesh.edof);

figure
for i=1:mesh.nel
   
    [stresses1, zcoords1, xCoord1, yCoord1] = el(i).computeStressThroughThickness(ed(:,i), [-1; 0]);
    [stresses2, zcoords2, xCoord2, yCoord2] = el(i).computeStressThroughThickness(ed(:,i), [1; 0]);
    [stressesMid, zcoordsMid, xCoordMid, yCoordMid] = el(i).computeStressThroughThickness(ed(:,i), [0; 0]);
    
    tmp = 1;
    plot_sigma = [stresses1.stress(1,tmp);    stresses2.stress(1,tmp)];
    plot_sigma_mid(i) = stressesMid.stress(1,tmp);
    plot_x_mid(i) = xCoordMid;
    
    plot_x = [xCoord1,xCoord2];
    plot(plot_x,plot_sigma,'k'); 
    hold on;
    drawnow;
end
% axis equal

ylabel('Stress [N/m^2]');
xlabel('x [m]');

plot([0 plot_x_mid], [-6e5 plot_sigma_mid])

% save('cb_stress_sse2','plot_sigma_mid','plotThetaMid')

