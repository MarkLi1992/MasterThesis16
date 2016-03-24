clear variables;

%Euler-Bernoulli vs SolidShell as h-> 0
P = -50; lx = 0.1, ly=0.01; E = 100e9; G = E/(2*(1+0));
h = linspace(1e-5,1e-2,100);

for ih = 1:length(h)
    nelx = 10; nely = 1; nlamel = 1;
    ss_maxdisp(ih) = solid8layered_runner('Konsol', 'SolidShell', h(ih), nelx, nely, nlamel);
    se_maxdisp(ih) = solid8layered_runner('Konsol', 'Solid', h(ih), nelx, nely, nlamel);

    Iy = (ly*h(ih)^3/12);
    ts_maxdisp(ih) = abs(P*lx/((5/6)*(ly*h(ih))*G) + P*lx^3/3/E/Iy);
    eb_maxdisp(ih) = abs(P*lx^3/3/E/Iy);
    aspectRatio(ih) = ly/h(ih);
end

plotXaxis = h;
figure
loglog(plotXaxis,eb_maxdisp);
hold on;
loglog(plotXaxis,ts_maxdisp,'--');
loglog(plotXaxis,ss_maxdisp);
loglog(plotXaxis,se_maxdisp);
legend('Euler-Bernoulli','Timoshenko','Solid shell', 'Solid')
xlabel('Thickness')
ylabel('displacement')