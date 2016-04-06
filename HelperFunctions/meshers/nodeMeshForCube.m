function [mesh] = nodeMeshForCube(nelx,nely,nelz, nnoxel,nnoyel,nnozel)              

nelno = nnoxel*nnoyel*nnozel;%Number of element nodes


totnnox = (nelx*(nnoxel-1)+1);
totnnoy = (nely*(nnoyel-1)+1);
totnnoz = (nelz*(nnozel-1)+1);

nno = totnnox*totnnoy*totnnoz;

nodelayout = reshape(1:nno, totnnox, totnnoy, totnnoz);

icel = 0;
for iz=0:(nnozel-1):totnnoz-nnozel
    for iy=0:(nnoyel-1):totnnoy-nnoyel
        for ix = 0:(nnoyel-1):totnnox-nnoyel
            icel = icel+1;
            temp = nodelayout((1:nnoxel) + ix,(1:nnoyel) +iy,(1:nnozel) + iz);
            temp = reshape(temp,1,nelno,1);
            mesh(icel,:) = temp;   
        end
    end
end
mesh = mesh';
