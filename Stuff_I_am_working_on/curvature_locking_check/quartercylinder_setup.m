function [ m, elprop, M, bc, ftrac ] = quartercylinder_setup( problem,nnn )

% Setup data for the problems
switch problem        
    case 'SR1'
         inner_radius=4.12; lx = 0.1; thickness = 0.2;
%         nelr = 2; nelphi=50; nelx=50;
%         nelr = 1; nelphi=60; nelx=6;
%         nelr = 5; nelphi=200; nelx=20;
%         nelr = 3; nelphi=100; nelx=10;
        nelr = nnn(1); nelphi=nnn(2); nelx=nnn(3);

        %Traction force
        tractionInfo.sides = [3];
        tractionInfo.tractions = {@(x,y) [0,0, -1/(thickness*lx)]'}; ap = (inner_radius*pi/nelphi)/(thickness/nelr);
        
        %angles
        angles = [0];
        
        %Material properties
%         EL = 174.6E9;
%         ET = 7E9;
%         nuLT = 0.25;
%         GLT = 3.5E9;
%         GTT = 1.4E9;
%         nuTL = ET/EL*nuLT;
%          
%         D = hooke_trans_iso(EL, ET, nuLT, GLT, GTT);
        D = hooke(4,10000000,0.25);
        
        m = Mesh();
        m.create_quarter_cylinder_mesh(lx, inner_radius, thickness, nelr, nelphi, nelx)

    otherwise
        error('Unknown problem type')
end

% Create element properties
elprop = ElementProperties();
elprop.setup(thickness, angles, D);

% Interpolation matrix
M = getInterPolMatrix(4);
% M = @(x,y,z)reshape([x,0.0,0.0,0.0,0.0,0.0,0.0,y,0.0,0.0,0.0,0.0,0.0,0.0,z,0.0,0.0,0.0,0.0,0.0,x.*z,0.0,0.0,0.0,0.0,0.0,y.*z,0.0,0.0,0.0,0.0,0.0,x.*y.*z,0.0,0.0,0.0,0.0,0.0,0.0,x,0.0,0.0,0.0,0.0,0.0,y,0.0,0.0,0.0,0.0,0.0,x.*y,0.0,0.0],[6,9]);
% Load and boundary conditions
bc = quartercylinderBC(problem, m.dof, m.sideElements);

% ftrac = zeros(m.ndofs,1);
% tmpdof = unique(m.dof(m.sideElements(3).nodes(:),3));
% ftrac(tmpdof) = -1/length(tmpdof);
% return;

ftrac = zeros(m.ndofs,1);
% ftrac = zeros(m.neldofs, m.nel);
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
            case 7
                e1 = ex; e2 = ey;
        end
        
        sideDofs = m.dof(cSideNodes(:,ii),1:3)';
        sideDofs = sideDofs(:);
        
        cElement = cSideElements(ii);
        %         feltrac = solid8trac(e1,e2,[0,0,-50/(lx*ly)]',1);
        feltrac = solid8traction(e1,e2,cTraction,3);
        
        ftrac(sideDofs) = ftrac(sideDofs) + feltrac;
%         ftrac(:,cElement) = ftemp(m.edof(:,cElement));
%         ftemp=zeros(m.ndofs,1);
    end
    
end



end

