function [bc] = cubeBC( opt, dof, sideElements)

%Get nodes for each side
side1 = sideElements(1).nodes;side1 = unique(side1(:));
side2 = sideElements(2).nodes;side2 = unique(side2(:));
side3 = sideElements(3).nodes;side3 = unique(side3(:));
side4 = sideElements(4).nodes;side4 = unique(side4(:));
side5 = sideElements(5).nodes;side5 = unique(side5(:));
side6 = sideElements(6).nodes;side6 = unique(side6(:));

if strcmp('Konsol', opt) || strcmp('Konsol_stacked',opt)
    ldofs = dof(side1,1:3);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
elseif strcmp('KonsolMedUtbredd', opt) || strcmp('KonsolMedUtbredd_stacked',opt)
    ldofs = dof(side1,:);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
elseif(strcmp('FrittUpplagt', opt))
    ldofs2 = dof(side2,2:3);
    ldofs2 = ldofs2(:);
    ldofs = dof(side1,:);

    ldofs = ldofs(:);
    bc = [[ldofs;ldofs2],[ldofs;ldofs2]*0];
elseif strcmp('InspandPlatta', opt) || strcmp('InspandPlatta_stacked',opt)
    ldofs = dof([side1;side2;side3;side4], 1:3);
    ldofs = ldofs(:);
    bc = [ldofs, ldofs*0];

elseif(strcmp('HybridStress2', opt))   ||   (strcmp('HybridStress2_stacked', opt))
    side1bottomDofs = dof(      sideElements(1).nodes(1:2,:)    ,1:3);
    side1bottomDofs = side1bottomDofs(:);
    side2bottomYZDofs = dof(      sideElements(2).nodes(1:2,:)    ,1:3);
    side2bottomYZDofs = side2bottomYZDofs(:);
    ldofs = [side1bottomDofs; side2bottomYZDofs];
%     ldofs = dof([side1],1:3); ldofs = ldofs(:);
%     ldofs = [ldofs; side2bottomYZDofs];
    bc = [ldofs, 0*ldofs];
   
elseif(strcmp('MembraneUnsymmetric', opt))
    ldofs = dof(side1,1:3);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
    
elseif(strcmp('CurvedBeam', opt)) || strcmp('CurvedBeam_stacked',opt)
    ldofs1 = dof(unique(side4),1:3);
    ldofs1 = ldofs1(:);
    bc = [ldofs1,ldofs1*0];
elseif(strcmp('NavierCheck', opt))
    side1 = sideElements(1).nodes(1:2,:); side1 = unique(side1(:));
    side2 = sideElements(2).nodes(1:2,:); side2 = unique(side2(:));
    side3 = sideElements(3).nodes(1:2,:); side3 = unique(side3(:));
    side4 = sideElements(4).nodes(1:2,:); side4 = unique(side4(:));
    
    bno = unique([side1;side2;side3;side4]);
    bdofsz = dof(bno,3);
    
    cornerDofs = dof(intersect(side1, side3),1:3); cornerDofs = cornerDofs(:);
    ldofs = [bdofsz; cornerDofs];
    
    bc = [ldofs,ldofs*0];
    elseif(strcmp('WhitneyCheck', opt))
    side1 = sideElements(1).nodes(1:2,:); side1 = unique(side1(:));
    side2 = sideElements(2).nodes(1:2,:); side2 = unique(side2(:));
    side3 = sideElements(3).nodes(1:2,:); side3 = unique(side3(:));
    side4 = sideElements(4).nodes(1:2,:); side4 = unique(side4(:));
    
    bno = unique([side1;side2;side3;side4]);
    bdofsz = dof(bno,3);
    
    s1x = dof(side1,1);
    s2x = dof(side2,1);
    s3y = dof(side3,2);
    s4y = dof(side4,2);
    
    ldofs = [bdofsz; s1x;s2x;s3y;s4y];
    
    bc = [ldofs,ldofs*0];
else
    error('nope');
end

end 

