function mesh = stressMesher(nLam, interp,ex,ey,ez, int_coordsG)
%
%Creates a mesh for a cube with 4 nodes in the plane, and X nodes in the
%out of plane direction, where X is defined for each stress component
%by the interp-variable.
%
%Input example:
% nLam = 2;
% interp = [2 2 3 2 4 4];

%Total number of dofs for each stress-component in the element
nno_s = (nLam + 1 + (interp-2)*nLam)*4;

%Number of stress componenets
ns = length(interp);

%number of dofs per layer
nno_s_lay = interp*4;

sdofs{1} = 1:nno_s(1);
for i = 2:ns
    sdofs{i} = (1:nno_s(i)) + max(sdofs{i-1});
end

for i=1:nLam    
    layDofs = [];
    for is=1:ns
        temp = sdofs{is};
        cDofs = temp((1:nno_s_lay(is)) + (nno_s_lay(is) - 4)*(i-1));
        layDofs = [layDofs;cDofs'];
    end
    edof(:,i) = layDofs;
end

mesh.edof = edof;
mesh.ndofs = max(max(edof));

mesh.ex = repmat(ex',1,nLam);
mesh.ey = repmat(ey',1,nLam);

for i=1:nLam
    mesh.ez(:,i) = [[1 1 1 1]*int_coordsG(i), [1 1 1 1]*int_coordsG(i+1)]';
end






