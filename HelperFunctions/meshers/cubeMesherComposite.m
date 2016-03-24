function [edof,coord,ex,ey,ez,dof,nel,ndofs,nno,...
    sideElements] = cubeMesherComposite(lx,ly,lz,nelx,nely,nlam,nlamel)
%Only works proporly for 2x2 nodes in-plane
% side1: back, 
% side2: front, 
% side3: left
% side4: right
% side5: top
%                         7
%                      ******
%                  ****   *  **
%               ****      *   **
%            ****         *     **
%         ***             *      **
%       8*                *       ***
%       ***               *         **
%       * **              *           *6
%       *   **            *          ***
%       *    **           *       ***  *
%       *      **         *    ***     *
%       *       **        *****        *
%       *         **    ***            *
%       *          **5**  *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    *            *
%       *            *    3            *
%       *            *  *****          *
%       *            ***    **         *
%       *        *****        **       *
%       *     ***    *         **      *
%       *  ***       *           **    *
%       ***          *            **   *
%       4*           *              ** *
%         **         *               ***
%          ***       *                *2
%            **      *             ***
%             **     *         ****
%               **   *      ****
%                **  *   ****
%                  ******
%                    1


ndofsno = 3;

nnoxel = 2;
nnoyel = 2;
nnozel = nlam+1 + (nlamel-1)*nlam;
nelz = 1;

nelno = nnoxel*nnoyel*nnozel;%Number of element nodes

elx = lx/(nelx);
ely = ly/nely;
elz = lz/nelz;

totnnox = (nelx*(nnoxel-1)+1);
totnnoy = (nely*(nnoyel-1)+1);
totnnoz = (nelz*(nnozel-1)+1);

xx = linspace(0,lx,totnnox);%0:elx:lx;
yy = linspace(0,ly,totnnoy);
zz = linspace(0,lz,totnnoz);

nno = totnnox*totnnoy*totnnoz;

ndofs = nno*ndofsno;
dof = reshape(1:ndofs,ndofsno,nno)';
nel = nelx*nely*nelz;

nodelayout = reshape(1:nno, totnnox, totnnoy, totnnoz);
ellayout = reshape(1:nel, nelx,nely,nelz);

icel = 0;
for iz=0:(nnozel-1):totnnoz-nnozel
    for iy=0:(nnoyel-1):totnnoy-nnoyel
        for ix = 0:(nnoyel-1):totnnox-nnoyel
            icel = icel+1;
            temp = nodelayout((1:nnoxel) + ix,(1:nnoyel) +iy,(1:nnozel) + iz);
            temp = reshape(temp,1,nelno,1); 
%             temp(3:4) = temp([4 3]); temp([7 8]) = temp([8 7]);
%             for itt=1:nnozel
%                 temp((3:4) + 4*(itt-1)) = temp(fliplr((3:4) + 4*(itt-1)))
%             end
            mesh(icel,:) = temp;
            edof(icel,:) = reshape(dof(temp,:)',1,nelno*ndofsno);
                
        end
    end
end
edof = edof';
mesh = mesh';

cn = 0;
for iz=1:totnnoz
    for iy=1:totnnoy
        for ix = 1:totnnox
            cn = cn +1;
            coord(cn,:) = [xx(ix), yy(iy), zz(iz)];
        end
    end
end
coord = coord';

xcoord = coord(1,:);
ycoord = coord(2,:);
zcoord = coord(3,:);

ex = xcoord(mesh);
ey = ycoord(mesh);
ez = zcoord(mesh);

if size(ex,1) == 1 & size(ex,2) > 1
   ex = ex'; ey = ey'; ez = ez'; 
end

side1nodes = nodelayout(1,:,:);
side2nodes = nodelayout(end,:,:);
side3nodes = nodelayout(:,1,:);
side4nodes = nodelayout(:,end,:);
side5nodes = nodelayout(:,:,end);
side6nodes = nodelayout(:,:,1);

side1el = ellayout(1,:,:);
side2el = ellayout(end,:,:);
side3el = ellayout(:,1,:);
side4el = ellayout(:,end,:);
side5el = ellayout(:,:,end);
side6el = ellayout(:,:,1);



%Side 1 elements
itr = 1;
% for iz = 1:nelz
for iz=0
    for iy = 1:nely
%         temp = side1nodes(1, (1:2) +iy-1,  (1:nnozel) +nnozel-1);
        temp = side1nodes(1, (1:2) +iy-1,  (1:nnozel));
        sideElements(1).elements(itr) = side1el(1,iy);
        sideElements(1).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end

%Side 2 elements
itr = 1;
% for iz = 1:nelz
for iz = 0
    for iy = 1:nely
%         temp = side2nodes(1, (1:2) +iy-1,  (1:2) +iz-1);
        temp = side2nodes(1, (1:2) +iy-1,  (1:nnozel));
%         sideElements(2).elements(itr) = side2el(iz,iy);
        sideElements(2).elements(itr) = side2el(1,iy);
        sideElements(2).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end

%Side 3 elements
itr = 1;
for iz = 0
    for ix = 1:nelx
        temp = side3nodes((1:2) +ix-1, 1, (1:nnozel));
        sideElements(3).elements(itr) = side3el(ix,1);
        sideElements(3).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end

%Side 4 elements
itr = 1;
for iz = 0
    for ix = 1:nelx
        temp = side4nodes((1:2) +ix-1, 1, (1:nnozel));
        sideElements(4).elements(itr) = side4el(ix,1);
        sideElements(4).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end

%Side 5 elements
itr = 1;
for ix = 1:nelx
    for iy = 1:nely
        temp = side5nodes((1:2) + ix-1, (1:2) +iy-1, 1);
        sideElements(5).elements(itr) = side5el(ix,iy);
        sideElements(5).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end

%Side 6 elements
itr = 1;
for ix = 1:nelx
    for iy = 1:nely
        temp = side6nodes((1:2) + ix-1, (1:2) +iy-1, 1);
        sideElements(6).elements(itr) = side6el(ix,iy);
        sideElements(6).nodes(:,itr) = temp(:);
        itr = itr+1;
    end
end


% outstruct.edof = edof;
% outstruct.coord = coord;
% outstruct.ex = ex;
% outstruct.ey = ey;
% outstruct.ez = ez;
% outstruct.dof = dof;
% outstruct.nel = nel;
% outstruct.ndofs = ndofs;
% outstruct.nno = nno;
% outstruct.ndofsel = nelno*3;
% outstruct.sideElements = sideElements;
% 
% outstruct.lx = lx;
% outstruct.ly = ly;
% outstruct.lz = lz;
% 
% outstruct.nelx = nelx;
% outstruct.nely = nely;
% outstruct.nelz = nelz;