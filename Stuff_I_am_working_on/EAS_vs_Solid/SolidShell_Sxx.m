clear variables

%Data for the "konsol-balk"
P = -50; lx = 0.1; ly=0.01; E = 100e9; G = E/(2*(1+0));
h = 0.002;
Iy = (ly*h^3/12);
% abs(P*lx^3/3/E/Iy)

% M_xx = @(x) P*lx*(1-x/lx);
q = P/(lx*ly)*ly; M_xx = @(x) q*lx*x - q*x.^2/2 - q*lx^2/2;

% sig_xx = @(x,z) (P*lx*(1-x/lx))  *z/Iy;
sig_xx = @(x,z) (M_xx(x))*z/Iy;

% Tforce = @(x) P;
Tforce = @(x) P/(lx*ly)*ly * (lx-x);

%Mesh options
nelx = 16; nely = 1; nlamel = 1;

aspectRatio = (lx/nelx)/(h/(nlamel))

%% Run
maxdisp = solid8layered_runner('KonsolMedUtbredd', 'SolidShell', h, nelx, nely, nlamel);

%% Postprocces
filename = sprintf('reference_konsol_sxx/SolidShell_x%i_y%i_z%i.mat',nelx,nely,nlamel);
load(filename);

%Get mid element
plotEl = coordinate2element(mesh.ex,mesh.ey,mesh.ez, [mesh.lx/2, mesh.ly/2, mesh.lz/2])

%DELETE______
iz =1;
for zz = linspace(0,h,40)
val(iz) = el(plotEl).computeTauIntegral(zz);
iz = iz+1;
end
figure;
plot(val, linspace(0,h,40));
%*******DELETE^^^^^^^^^^


%Plot shear force
figure
plot([0,lx], Tforce([0,lx]));
hold on;
for plotEl = 1:nelx
[T1, outcoords1] = el(plotEl).computeShearForceAtXY(a(mesh.edof(:,plotEl)),[-1,0], [0, h]);
[T2, outcoords2] = el(plotEl).computeShearForceAtXY(a(mesh.edof(:,plotEl)),[1,0], [0, h]);
T1 = T1*ly;
T2 = T2*ly;
plot([outcoords1(1),outcoords2(1)], [T1,T2],'r*-');
end
% axis equal

%Plot moments
figure
plot([0,lx], M_xx([0,lx]));
hold on;
for plotEl = 1:nelx
[M1, outcoords1] = el(plotEl).computeMomentAtXY(a(mesh.edof(:,plotEl)),[-1,0], [0, h]);
[M2, outcoords2] = el(plotEl).computeMomentAtXY(a(mesh.edof(:,plotEl)),[1,0], [0, h]);
M1 = M1*ly;
M2 = M2*ly;
plot([outcoords1(1),outcoords2(1)], [M1,M2],'r*-');
end

%Plot stresses
stressesAtPoints = [-(1/sqrt(3)), (1/sqrt(3));...
    0             0;...
    -1           -1];
figure
analytXX = linspace(0,lx,100);
plot(analytXX, sig_xx(analytXX,-h/2));
hold on; xlabel('x-direction'); ylabel('\sigma_{xx}');
for plotEl = 1:nelx
    [s1, gpcoords1] = el(plotEl).computeStressAt(a(mesh.edof(:,plotEl)),[-1 0 -1]');
    [s2, gpcoords2] = el(plotEl).computeStressAt(a(mesh.edof(:,plotEl)),[1 0 -1]');
    
    %% Plot
    plot([gpcoords1(1), gpcoords2(1)], [s1(1).stress(1), s2(1).stress(1)],'r-*');
end




