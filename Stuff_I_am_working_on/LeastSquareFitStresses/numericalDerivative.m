function [deriv] = numericalDerivative(funk, x0)
% h = 1e-7;
% deriv = (funk(x0+h) - funk(x0))/h;
nvars = length(x0);
temp = zeros(nvars,1);
opt = 1;
if(opt == 1)
    
    h = 1e-20;
    for ii=1:nvars
        pert = temp;
        pert(ii) = h;
        
        deriv(:,ii) = imag( funk(x0 + 1i*pert)  )/h;
    end
    
else
    
    h = 1e-8;
    f0 = funk(x0);
    for ii=1:nvars
        pert = temp;
        pert(ii) = h;
        
        deriv(:,ii) = ( funk(x0 + pert) - f0 )/h;
    end
end

end



