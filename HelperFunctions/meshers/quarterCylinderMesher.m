function [edof,coord,ex,ey,ez,dof,nel,ndofs,nno,...
    mesh, sideElements] = quarterCylinderMesher( length, inner_radius, thickness, nelr, nelphi, nelx, el_nnor, el_nnophi, el_nnox, ndofsno)

%Mesh conectivity is the same as for a plate, but with fused side nodes
[edof,coord,~,~,~,dof,nel,ndofs,nno,...
    mesh, sideElements] = cubeMesher(1,1,1, nelx,nelphi,nelr, el_nnox,el_nnophi, el_nnor, ndofsno);

outer_radius = inner_radius + thickness;

totnnophi = (nelphi*(el_nnophi-1)+1);
totnnox = (nelx*(el_nnox-1)+1);
totnnor = (nelr*(el_nnor-1)+1);

cx = linspace(0,length,totnnox);
cr = linspace(inner_radius,outer_radius,totnnor);
ang_dif = (pi/2)/(totnnophi-1);
currendnode = 1; 

for ir = 1:totnnor
    currentRadius = cr(ir);
    acc_ang = -ang_dif;
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

xcoord = coord(1,:);
ycoord = coord(2,:);
zcoord = coord(3,:);

ex = xcoord(mesh);
ey = ycoord(mesh);
ez = zcoord(mesh);

st1 = ceil(nelphi/2);
en1 = nelphi*(nelx - 1) + st1;

sideElements(7).elements = sideElements(5).elements(st1:nelphi:en1 );
sideElements(7).nodes    = sideElements(5).nodes(:, st1:nelphi:en1);

end

