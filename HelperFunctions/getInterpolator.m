function [ stressInterp ] = getInterpolator( interpx,interpy,interpz )


for ii=1:length(interpz)
    ix = num2str(interpx(ii));
    iy = num2str(interpy(ii));
    iz = num2str(interpz(ii));
    eval(['stressInterp{ii} = InterpolatorX',ix,'Y',iy,'Z',iz,';']); 
end

end

