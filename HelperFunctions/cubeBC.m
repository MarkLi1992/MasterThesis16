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
elseif(strcmp('InspandPlatta', opt))
    ldofs = dof([side1;side2;side3;side4], 1:3);
    ldofs = ldofs(:);
    bc = [ldofs, ldofs*0];

elseif(strcmp('HybridStress2', opt))
    side1bottomDofs = dof(      sideElements(1).nodes(1,:)    ,1:3);
    side1bottomDofs = side1bottomDofs(:);
    side2bottomYZDofs = dof(      sideElements(2).nodes(1,:)    ,2:3);
    side2bottomYZDofs = side2bottomYZDofs(:);
    bc = [[side1bottomDofs; side2bottomYZDofs], 0*[side1bottomDofs; side2bottomYZDofs]];
elseif(strcmp('Testdelete', opt))
    dofsu = dof( intersect(side1,side5),1:3); dofsu = dofsu(:);
    dofsl = dof( intersect(side1,side6),1:3); dofsl = dofsl(:);
%     bc = [dofsu;dofsl];
%     bc = [bc,bc*0];
    
    dofsl2 = dof( intersect(side2,side6),2:3); dofsl2 = dofsl2(:);
    
    bc = [dofsl;dofsl2];
    bc = [bc,bc*0];
    

elseif(strcmp('MembraneUnsymmetric', opt))
    ldofs = dof(side1,1:3);
    ldofs = ldofs(:);
    bc = [ldofs,ldofs*0];
    
elseif(strcmp('TestUdiscont', opt))
    global stepp
    if(stepp == 1)
        ldofs = dof(side1,:);
        ldofs = ldofs(:);
        bc = [ldofs,ldofs*0];
    else
        global eldisp
        lowdofs = dof(side6,:)';
        updofs  = dof(side5,:)';
        dofs = [lowdofs(:);updofs(:)];
        bc = [dofs, eldisp];
    end
else
    error('nope');
end

end 

