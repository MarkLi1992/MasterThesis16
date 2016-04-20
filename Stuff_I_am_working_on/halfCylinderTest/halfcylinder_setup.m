function [ m, elprop, M, bc, ftrac ] = halfcylinder_setup( problem,nnn )

% Setup data for the problems
switch problem        
    case 'SR1'
         inner_radius=0.04; lx = 0.03*pi; thickness = 0.002;
%         nelr = 2; nelphi=50; nelx=50;
%         nelr = 1; nelphi=50; nelx=50;
        nelr = nnn(1); nelphi=nnn(2); nelx=nnn(3);

        %Traction force
        tractionInfo.sides = [7];
        tractionInfo.tractions = {@(x,y) [0,0, -500000]'}; ap = (lx/nelx)/(thickness/nelr);
        
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
        D = hooke(4,100e9,0.3);
        
        m = Mesh();
        m.create_half_cylinder_mesh(lx, inner_radius, thickness, nelr, nelphi, nelx)

    otherwise
        error('Unknown problem type')
end

% Create element properties
elprop = ElementProperties();
elprop.setup(thickness, angles, D);

% Interpolation matrix
M = getInterPolMatrix(6);

% Load and boundary conditions
bc = halvcylinderBC(problem, m.dof, m.sideElements);

ftrac = zeros(m.ndofs,1);
tmpdof = unique(m.dof(m.sideElements(7).nodes(:),3));
ftrac(tmpdof) = -5000/length(tmpdof);
return;
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

