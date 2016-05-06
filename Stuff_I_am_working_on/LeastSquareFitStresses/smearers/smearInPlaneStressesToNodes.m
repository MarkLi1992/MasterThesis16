function [as] = smearInPlaneStressesToNodes(mesh, el, a)
% components = [1 2 3 4 5 6];   
components = [1 2 4];
ncmp = length(components);

ed = a(mesh.edof);
sig = zeros(mesh.nno,ncmp,el(1).elprop.nLam);
counter = zeros(mesh.nno,1);
for iel = 1:mesh.nel
    elementNodes = mesh.nomesh(:,iel);
    [stresses] = el(iel).computeStressThroughLayers(ed(:,iel), [0,0]);
    
    for ilay = 1:size(stresses,3)
        %For each componenet xx, yy, xy
%         components = [1 2 4];
%         components = [1 2 3 4 5 6];
        for ic = 1:ncmp
            sig(elementNodes(1:4),ic,ilay) = sig(elementNodes(1:4),ic,ilay) + stresses(components(ic),1,ilay); %row = sigma_ic, col = bottom or top, i:th layer
            sig(elementNodes(5:8),ic,ilay) = sig(elementNodes(5:8),ic,ilay) + stresses(components(ic),2,ilay);
        end
    end
    counter(elementNodes) = counter(elementNodes) + 1;

end

as = zeros(mesh.nno*ncmp,el(1).elprop.nLam);

for ii = 1:el(1).elprop.nLam
   
     smeared = sig(:,:,ii) ./  repmat(counter,1,ncmp);
     for ic = 1:ncmp
         as(ic:ncmp:end,ii) = smeared(:,ic);
     end
%      as(1:3:end,ii) = smeared(:,1);
%      as(2:3:end,ii) = smeared(:,2);
%      as(3:3:end,ii) = smeared(:,3);
end


