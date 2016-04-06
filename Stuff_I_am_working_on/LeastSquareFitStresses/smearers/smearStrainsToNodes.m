function [as] = smearStrainsToNodes(mesh, el, a)

ed = a(mesh.edof);
nodes.epsxx = zeros(mesh.nno,1);
nodes.epsyy = zeros(mesh.nno,1);
nodes.epsxy = zeros(mesh.nno,1);
nodes.counter = zeros(mesh.nno,1);
for iel = 1:mesh.nel
   
    elementNodes = mesh.nomesh(:,iel);
    [~,stop,~] = el(iel).computeStressAt(ed(:,iel), [0,0,1]);
    [~,sbot,~] = el(iel).computeStressAt(ed(:,iel), [0,0,-1]);
    
    nodes.epsxx(elementNodes(1:4)) = nodes.epsxx(elementNodes(1:4)) + sbot(1);
    nodes.epsxx(elementNodes(5:8)) = nodes.epsxx(elementNodes(5:8)) + stop(1);
    
    nodes.epsyy(elementNodes(1:4)) = nodes.epsyy(elementNodes(1:4)) + sbot(2);
    nodes.epsyy(elementNodes(5:8)) = nodes.epsyy(elementNodes(5:8)) + stop(2);
    
    nodes.epsxy(elementNodes(1:4)) = nodes.epsxy(elementNodes(1:4)) + sbot(4);
    nodes.epsxy(elementNodes(5:8)) = nodes.epsxy(elementNodes(5:8)) + stop(4);
    
    nodes.counter(elementNodes) = nodes.counter(elementNodes) + 1;
    
end

as = zeros(mesh.nno*3,1);

sxx = nodes.epsxx./nodes.counter;
syy = nodes.epsyy./nodes.counter;
sxy = nodes.epsxy./nodes.counter;

as(1:3:end) = sxx;
as(2:3:end) = syy;
as(3:3:end) = sxy;


