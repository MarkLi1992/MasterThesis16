function [m, elprop, M, bc, ftrac ] = setup_problem(problem)
%setup problem: creates all the necessary data for a given problem

% Corresponds to an 8 noded brick element with 3 dofs per node


% Setup data for the problems
switch problem
    
    case 'InspandPlatta'
        % Mesh parameters
        lx = 0.6; ly=0.6; lz = 0.005;
        nelx = 40; nely=40; nelz=1;
        
        % Loading
        tractionInfo.sides = [5];
        tractionInfo.tractions = {@(x,y) [0 0 -5e3]'};
        
        %Material properties
        EL = 50E9;    % [Pa]
        ET = 9E9;   % [Pa]
        nuLT = 0.22;    % [-]
        GLT = 5E9;  % [Pa]
        GTT = 3.2E9;     % [Pa]
        
        D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT );
        %         D = hooke(4, 100e9, 0);
        angles = [0 90 0];
        
    case 'KonsolMedUtbredd'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 10; nely=1; nelz=1;
        
        %Traction force
        tractionInfo.sides = [5];
        tractionInfo.tractions = {@(x,y) [0 0 -50/(lx*ly)]'};
        
        %angles
        angles = [0 90 0 90];
        
        %Material properties
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
        D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
%         D = hooke(4,100e9,0);
    case 'HybridStress2'
        lx = 0.1; ly=0.01; lz = 0.01;
        nelx = 10; nely=1; nelz=1;

        tractionInfo.sides = [5];
        tractionInfo.tractions = {@(x,y) [0 0 -5/(lx*ly)]'};%*sin(pi*y/ly)
        
        angles = [0 90];
        %Material properties
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
%         D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
D = hooke(4,100e9,0);
    case 'Konsol'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 20; nely=1; nelz=1;
        
        tractionInfo.sides = [2];
        tractionInfo.tractions = {@(x,y) [0,0, -50/(lz*ly)]'};
        
        angles = [0];
        
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
%         D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
                D = hooke(4,100e9,0);
        
    case 'MembraneUnsymmetric'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 10; nely=2; nelz=1;
        
        tractionInfo.sides = [2];
        tractionInfo.tractions = {@(x,y) [1000 0 0]'};
        
        angles = [0 0 90 90];
        
        %Material properties
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
        D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT );
        
    case 'Konsol_stacked'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 20; nely=1; nlamel=1;
        
        tractionInfo.sides = [2];
        tractionInfo.tractions = {@(x,y) [0 0 -50/(lz*ly)]'};
        
        angles = [0];
        nlam = length(angles);
        D = hooke(4,100e9,0);
        
        m = Mesh();
        m.create_cube_mesh_stacked_solid_elements(lx,ly,lz,nelx,nely,nlam, nlamel);
        
    case 'KonsolMedUtbredd_stacked'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 20; nely=2; nlamel=6;
        
        tractionInfo.sides = [5];
        tractionInfo.tractions = {@(x,y) [0 0 -50/(lx*ly)]'};
        
        angles = [0];
        nlam = length(angles);
        %Material properties
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
        D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
        
        %independent mesh
        m = Mesh();
        m.create_cube_mesh_stacked_solid_elements(lx,ly,lz,nelx,nely,nlam, nlamel);
        
    otherwise
        error('Unknown problem type')
end

% Create element properties
elprop = ElementProperties();
elprop.setup(lz, angles, D);

% Create mesh
if ~(exist('m','var'))
m = Mesh();
m.create_cube_mesh(lx,ly,lz,nelx,nely,nelz);
end

% Interpolation matrix
M = getInterPolMatrix(1);

% Load and boundary conditions
bc = cubeBC(problem, m.dof, m.sideElements);

ftemp = zeros(m.ndofs,1);
ftrac = zeros(m.neldofs, m.nel);
%setup traction-force-vector
for iside = 1:length(tractionInfo)
    cSideIndex = tractionInfo.sides(iside);
    cSideNodes = m.sideElements(cSideIndex).nodes;
    cSideElements = m.sideElements(cSideIndex).elements;
    
    cTraction =  tractionInfo.tractions{iside};
    
    for ii = 1:size(cSideNodes,2)
        
        sideNodes = cSideNodes(:,ii);
        
        %Simple hack to get to only get for nodes (Needed if
        %sidesWithTraction is not equal to 5 or 6)
        sideNodes = sideNodes([1 2 end-1, end]);
        
        ncoords = m.coord(:,sideNodes);
        ex = ncoords(1,:); ey = ncoords(2,:); ez = ncoords(3,:);
        
        switch cSideIndex
            case 1
                e1 = ey; e2 = ez;
            case 2
                e1 = ey; e2 = ez;
            case 3
                e1 = ex; e2 = ez;
            case 4
                e1 = ex; e2 = ez;
            case 5
                e1 = ex; e2 = ey;
                case 6
                e1 = ex; e2 = ey;
        end
        
        sideDofs = m.dof(cSideNodes(:,ii),1:3)';
        sideDofs = sideDofs(:);
        
        cElement = cSideElements(ii);
        %         feltrac = solid8trac(e1,e2,[0,0,-50/(lx*ly)]',1);
        feltrac = solid8traction(e1,e2,cTraction,3);
        
        ftemp(sideDofs) = feltrac;
        ftrac(:,cElement) = ftemp(m.edof(:,cElement));
        ftemp=zeros(m.ndofs,1);
    end
    
end

