function [m, elprop, M, bc, ftrac, nlamel] = setup_problem_solidComposite(problem, thickness, nelx, nely, nlamel)
%setup problem: creates all the necessary data for a given problem

% Corresponds to an 8 noded brick element with 3 dofs per node


% Setup data for the problems
switch problem
    
    case 'InspandPlatta'
        % Mesh parameters
        lx = 0.6; ly=0.3; lz = thickness;
%         nelx = 50; nely=50;
%         nlamel = 5; %Number of elements per layer

        % Loadings
        sidesWithTraction = 5;
        sideTractions = [0 0 -5000]';
        
        %Material properties
        EL = 50E9;    % [Pa]
        ET = 9E9;   % [Pa]
        nuLT = 0.22;    % [-]
        GLT = 5E9;  % [Pa]
        GTT = 3.2E9;     % [Pa]
        
        D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT );
        angles = [0 90 0];
        
    case 'Konsol'
        lx = 0.1; ly=0.01; lz = thickness;
%         nelx = 10; nely=1; 
%         nlamel = 1; %Number of elements per layer
        
        sidesWithTraction = 2;
        sideTractions = [0 0 -50/(lz*ly)]';

        angles = [0];
        
%         EL = 174.6E9;
%         ET = 7E9;
%         nuLT = 0.25;
%         GLT = 3.5E9;
%         GTT = 1.4E9;
%         nuTL = ET/EL*nuLT;
        
%         D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
        D = hooke(4,100e9,0);

        aspecRatio = (lx/nelx)/(lz/(nlamel*length(angles)));
    case 'KonsolMedUtbredd'
        lx = 0.1; ly=0.01; lz = thickness;
%         nelx = 10; nely=1; 
%         nlamel = 1; %Number of elements per layer
        
        
        sidesWithTraction = 5;
        sideTractions = [0 0 -50/(lx*ly)]';

        angles = [0];
        aspecRatio = (lx/nelx)/(lz/(nlamel*length(angles)));
        
        D = hooke(4,100e9,0);
        %Material properties
%         EL = 50E9;
%         ET = 9E9;
%         nuLT = 0.22;
%         GLT = 5E9;
%         GTT = 3.2E9;
%         nuTL = ET/EL*nuLT;
%         
%         D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT );
%         
    case 'HybridStress2'
        %A         B       H
        lx = 0.1; ly=0.01; lz = thickness;
%         nelx = 10; nely=1; 
%         nlamel = 1; %Number of elements per layer

%         lx/lz
        sidesWithTraction = [5];
        sideTractions = [0 0 -50/(lx*ly)]';%*sin(pi*y/ly)

        angles = [0 90];
        %Material properties
        EL = 174.6E9;
        ET = 7E9;
        nuLT = 0.25;
        GLT = 3.5E9;
        GTT = 1.4E9;
        nuTL = ET/EL*nuLT;
        
        D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
    case 'TestUdiscont'
        global stepp;
        if(stepp == 1)
            lx = 0.1; ly=0.01; lz = 0.005;
            nelx = 10; nely=1; 
            nlamel = 1; %Number of elements per layer
            
                    sidesWithTraction = 5;
                sideTractions = [0 0 -50/(lx*ly)]';
        else
            lx = 0.1/10; ly=0.01; lz = 0.005;
            nelx = 1; nely=1; 
            nlamel = 2; %Number of elements per layer
            
                    sidesWithTraction = 5;
                sideTractions = [0 0 -0/(lx*ly)]';
        end
        


        angles = [0];
        aspecRatio = (lx/nelx)/(lz/(nlamel*length(angles)));
        
        D = hooke(4,100e9,0);

    otherwise
        error('Unknown problem type')        
end

% Create element properties
elprop = ElementProperties();
elprop.setup(lz, angles, D);

% Create mesh
m = Mesh;
m.create_cube_mesh_stacked_solid_elements(lx,ly,lz,nelx,nely,elprop.nLam,nlamel)
% m = cubeMesherComposite(lx,ly,lz,nelx,nely,elprop.nLam,nlamel);

% Interpolation matrix
M = getInterPolMatrix(1);

% Load and boundary conditions
bc = cubeBC(problem, m.dof, m.sideElements);
ftrac=zeros(m.ndofs,1);

%setup traction-force-vector
for iside = 1:length(sidesWithTraction)
    currentSideIndex = sidesWithTraction(iside);
    currentSideNodes = m.sideElements(currentSideIndex).nodes;
    currentSideElements = m.sideElements(currentSideIndex).elements;
    currentTraction =  sideTractions(:,iside);

    for ii = 1:size(currentSideNodes,2)
        
        sideNodes = currentSideNodes(:,ii);
        
        %Simple hack to get to only get for nodes (Needed if
        %sidesWithTraction is not equal to 5 or 6)
        sideNodes = sideNodes([1 2 end-1, end]);
        
        %Get element-coords
        ncoords = m.coord(:,sideNodes);
        ex = ncoords(1,:); ey = ncoords(2,:); ez = ncoords(3,:);
        
        switch currentSideIndex
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
        end
        t = currentTraction ;

        sideDofs = m.dof(sideNodes,1:3)';
        sideDofs = sideDofs(:);

        feltrac = solid8trac(e1,e2,t,1);

        ftrac(sideDofs) = ftrac(sideDofs) + feltrac;
    end
    
end

