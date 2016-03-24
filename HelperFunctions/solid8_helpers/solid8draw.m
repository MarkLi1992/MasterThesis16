function [ ] = solid8draw( ex,ey,ez, nnodes)

if (~exist('nnodes','var'))
   nnodes = 8; 
end

nel = size(ex,2);
hold on;
for ie = 1:nel
    
   for it = 1:nnodes 
     cel(it,1:3) = [ex(it,ie), ey(it,ie), ez(it,ie)];  
   end
   
   if nnodes == 8
       %    order = [1 2 3 4 1 5 8 7 6 5 8 4 3 7 6 2];
      order = [1 2 4 3 1 5 7 8 6 5 7 3 4 8 6 2];
   elseif nnodes == 27
       order = [1 2 3 6 9 8 7 4 1 [1 2 3 6 9 8 7 4 1]+9 [1 2 3 6 9 8 7 4 1]+18];
   else
       error('nnodes error')
   end

   
   drawx = cel(order,1);
   drawy = cel(order,2);
   drawz = cel(order,3);
   
   plot3(drawx,drawy,drawz,'k'); 
   
end
hold off;
end

