function [ out ] = smearToNodes_MultipleLayers( nomesh, es, nodesInThickness )

nlay = size(es,3);

out = [];
for i=1:nlay
   as = smearToNodes(nomesh, es(:,:,i), nodesInThickness);
   out = [out, as];
end


end

