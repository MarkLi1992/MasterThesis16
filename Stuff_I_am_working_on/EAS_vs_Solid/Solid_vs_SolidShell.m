clear variables;

%Data for the "konsol-balk"
P = -50; lx = 0.1; ly=0.01; E = 100e9; G = E/(2*(1+0));
% h = 0.005;
h = 0.01;
Iy = (ly*h^3/12);
ts_maxdisp = abs(P*lx/((5/6)*(ly*h)*G) + P*lx^3/3/E/Iy);
eb_maxdisp = abs(P*lx^3/3/E/Iy);
tauxz = @(x,z,b,h,L,q) 3/2/(b*h)*(1 - (z/(h/2)).^2)*q;

%Mesh options 
nelx = 10; nely = 1; nlamel = 10;
% nelx = 10; nely = 1; nlamel = 10;

aspectRatio = (lx/nelx)/(h/(nlamel))

%Run that script
% maxdisp = solid8layered_runner('Konsol', 'SolidShell', h, nelx, nely, nlamel)
maxdisp = solid8layered_runner('Konsol', 'Solid', h, nelx, nely, nlamel)

%% Postprocces that script
% filename = 'reference_solution_konsol/Solid_x201_y1_z10.mat';
filename = sprintf('reference_solution_konsol/Solid_x%i_y%i_z%i.mat',nelx,nely,nlamel);
% filename = 'reference_solution_konsol/Solid_x100_y1_z10.mat';
[stresses, zcoords, xCoord] = solid8layered_postprocess(filename);

%% Plot da shitz
stressesAtPoints = [-1, 1 1 -1 0; -1 -1 1 1 0]*(1/sqrt(3));

%Txz
stressComp = 5;
figure; title('\Tau_{xz}')
for ip=1:size(stressesAtPoints,2)
    analytZZ = linspace(-h/2,h/2,20);
    analytTauxz = tauxz(xCoord(ip),analytZZ,ly,h,lx,-50);
    subplot(2,3,ip)
    plot(stresses(ip).stress(stressComp,:), zcoords);
    hold on;
    plot(analytTauxz, analytZZ + h/2);
    title(sprintf('x=%.3f, y=%.3f', stressesAtPoints(1,ip), stressesAtPoints(2,ip) ))
end












