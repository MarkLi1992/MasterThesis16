function [as] = smearToNodes(nomesh, es, nodesInThickness)

nel = size(es,2);
nno = max(max(nomesh));
totalValues = zeros(nno,1);
counter = zeros(nno,1);
for iel = 1:nel
    celementNodes = nomesh(:,iel);
    cvalues = es(:,iel);
    
    for i=1:nodesInThickness
       temp = (1:4) + ((i-1)*4);
       totalValues(celementNodes(temp)) = totalValues(celementNodes(temp)) + cvalues(i);
    end
    
    counter(celementNodes) = counter(celementNodes) + 1;

end

as = totalValues./counter;


