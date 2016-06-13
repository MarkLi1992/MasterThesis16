
ed = a(mesh.edof);

figure
for i=1:mesh.nel
   
    [stresses1, zcoords1, xCoord1, yCoord1] = el(i).computeStressThroughThickness(ed(:,i), [0; -1]);
    [stresses2, zcoords2, xCoord2, yCoord2] = el(i).computeStressThroughThickness(ed(:,i), [0; 1]);
    [stressesMid, zcoordsMid, xCoordMid, yCoordMid] = el(i).computeStressThroughThickness(ed(:,i), [0; 0]);
    
    tmp = 1;
    plot_sigma = [stresses1.stress(2,tmp);    stresses2.stress(2,tmp)];
    plot_sigma_mid(i) = stressesMid.stress(2,tmp);
   
    plot_x = [i, i+1];
    radius = sqrt(   yCoordMid^2 + mean(zcoordsMid(1))^2  );
    plotThetaMid(i) = asin( yCoordMid / radius);
    
%     plot(plot_x,plot_sigma);
%     hold on;
    plot(plotThetaMid(i),plot_sigma_mid(i),'*'); hold on;
    drawnow;
end


% save('cb_stress_sse2','plot_sigma_mid','plotThetaMid')

