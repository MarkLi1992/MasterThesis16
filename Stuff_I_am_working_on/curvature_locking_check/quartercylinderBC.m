function [bc] = quartercylinderBC( opt, dof, sideElements)

%Get nodes for each side
side1 = sideElements(1).nodes;side1 = unique(side1(:));
side2 = sideElements(2).nodes;side2 = unique(side2(:));
side3 = sideElements(3).nodes;side3 = unique(side3(:));
side4 = sideElements(4).nodes;side4 = unique(side4(:));
side5 = sideElements(5).nodes;side5 = unique(side5(:));
side6 = sideElements(6).nodes;side6 = unique(side6(:));

switch opt
    case 'SR1'
        ldofs1 = dof(unique(side4),1:3);
        ldofs1 = ldofs1(:);
        bc = [ldofs1,ldofs1*0];
    otherwise
    error('nope');
end

end 

