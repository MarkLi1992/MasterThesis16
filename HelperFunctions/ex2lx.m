function [lx,ly,lz] = ex2lx(ex,ey,ez, layerz)

nlay = length(layerz)-1;
laythick = layerz(end);

ncorners = 8;
connectedCorners = [1 5; 2 6; 3,7; 4,8];
cornCoord = zeros(ncorners,3);
lz = zeros(8,nlay); lx = lz; ly = lz;

for i=1:size(connectedCorners,1)
   ci1 = connectedCorners(i,1);
   ci2 = connectedCorners(i,2);
   cornCoord(ci1,:) = [ex(ci1),ey(ci1),ez(ci1)];
   cornCoord(ci2,:) = [ex(ci2),ey(ci2),ez(ci2)];
    
   thick = norm(cornCoord(ci2,:)-cornCoord(ci1,:));
   scl = thick/laythick;
   layerz_e = layerz*scl;
   
   dir = (cornCoord(ci2,:)-cornCoord(ci1,:))/thick;
   
   for jj = 1:nlay
       tmp1 = dir*layerz_e(jj);
       tmp2 = dir*layerz_e(jj+1);
       lx([ci1,ci2],jj) = cornCoord(ci1,1) + [tmp1(1),tmp2(1)];
       ly([ci1,ci2],jj) = cornCoord(ci1,2) + [tmp1(2),tmp2(2)];
       lz([ci1,ci2],jj) = cornCoord(ci1,3) + [tmp1(3),tmp2(3)];
   end
   
   %save variables for check
   savethick(i) = thick;
end


% keyboard;

end