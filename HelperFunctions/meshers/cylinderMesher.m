function [edof,coord,ex,ey,ez,dof,nel,ndofs,nno,...
    mesh, sideElements] = cylinderMesher( length, inner_radius, thickness, nelr, nelphi, nelx, el_nnor, el_nnophi, el_nnox, ndofsno)

%Mesh conectivity is the same as for a plate, but with fused side nodes
[~,~,~,~,~,~,~,~,~,...
    mesh, sideElementsCube] = cubeMesher(1,1,1, nelx,nelphi,nelr, el_nnox,el_nnophi, el_nnor, ndofsno);

outer_radius = inner_radius + thickness;

totnnophi = (nelphi*(el_nnophi-1)+1) - 1;
totnnox = (nelx*(el_nnox-1)+1);
totnnor = (nelr*(el_nnor-1)+1);
nelno = el_nnor*el_nnophi*el_nnox;

nno = totnnophi*totnnox*totnnor;
ndofs = nno * ndofsno;
dof = reshape(1:ndofs,ndofsno,nno)';
nel = nelr*nelx*nelphi;

cx = linspace(0,length,totnnox);
cr = linspace(inner_radius,outer_radius,totnnor);
ang_dif = (2*pi)/totnnophi;
currendnode = 1; 

for ir = 1:totnnor
    currentRadius = cr(ir);
    acc_ang = (-pi/4);
    for i=1:totnnophi
        acc_ang = acc_ang + ang_dif;
        y = cos(acc_ang)*currentRadius;
        z = sin(acc_ang)*currentRadius;
        
        for ix=1:totnnox
            coord(:,currendnode) = [cx(ix),z,y]';
            currendnode = currendnode +1;
        end
    end
end

side4nodes = sideElementsCube(4).nodes(:);
side3nodes = sideElementsCube(3).nodes(:);
for i = 1:size(side4nodes)
    mesh(side4nodes(i) == mesh) = side3nodes(i); 

    sideElementsCube(1).nodes(sideElementsCube(1).nodes == side4nodes(i)) = side3nodes(i); 
    sideElementsCube(2).nodes(sideElementsCube(2).nodes == side4nodes(i)) = side3nodes(i); 
    sideElementsCube(5).nodes(sideElementsCube(5).nodes == side4nodes(i)) = side3nodes(i); 
    sideElementsCube(6).nodes(sideElementsCube(6).nodes == side4nodes(i)) = side3nodes(i); 
end

count = 1; inum = 1;
while inum <=nno
    if   sum(sum(mesh==count)) > 0 
        mesh(mesh == count) = inum;
        sideElementsCube(1).nodes(sideElementsCube(1).nodes == count)  = inum;
        sideElementsCube(2).nodes(sideElementsCube(2).nodes == count)  = inum;
        sideElementsCube(5).nodes(sideElementsCube(5).nodes == count)  = inum;
        sideElementsCube(6).nodes(sideElementsCube(6).nodes == count)  = inum;
        inum = inum + 1;
    end
    count = count +1;
end

sideElements(1).nodes = sideElementsCube(1).nodes;
sideElements(1).elements = sideElementsCube(1).elements;
sideElements(2).nodes = sideElementsCube(2).nodes;
sideElements(2).elements = sideElementsCube(2).elements;
sideElements(3).nodes = sideElementsCube(6).nodes;
sideElements(3).elements = sideElementsCube(6).elements;
sideElements(4).nodes = sideElementsCube(5).nodes;
sideElements(4).elements = sideElementsCube(5).elements;


for i=1:nel
    temp = mesh(:,i);
    edof(:,i) = reshape(dof(temp,:)',1,nelno*ndofsno);
end

xcoord = coord(1,:);
ycoord = coord(2,:);
zcoord = coord(3,:);

ex = xcoord(mesh);
ey = ycoord(mesh);
ez = zcoord(mesh);

end

