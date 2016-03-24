function el = coordinate2element(ex,ey,ez,coord)

nel = size(ex,2);
for ii = 1:nel
    if  (max(ex(:,ii))>=coord(1) &&  min(ex(:,ii))<= coord(1)) &&...
            (max(ey(:,ii))>=coord(2) &&  min(ey(:,ii))<= coord(2)) &&...
            (max(ez(:,ii))>=coord(3) &&  min(ez(:,ii))<= coord(3))
        
        el = ii;
        return;
    end
    
end

end

