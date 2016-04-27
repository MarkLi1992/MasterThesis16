function [M] = getInterPolMatrix(opt)

if(opt == 1)

    M = @(xi,eta,zeta) [[0 0 0 0]; [0 0 0 0]; [0 0 0 0] ;zeros(1,4); [xi,0 xi*eta, 0; 0,eta,0,xi*eta]];
%     M = @(xi,eta,zeta) [[0 0 0 0 0]; [0 0 0 0 0]; [0 0 0 0 zeta] ;zeros(1,5); [xi,0 xi*eta, 0 0; 0,eta,0,xi*eta 0]];
%     M = @(xi,eta,zeta) [[0 0 0 0 xi];zeros(3,5); [xi,0 xi*eta, 0 0; 0,eta,0,xi*eta 0]];
    
%     M = @(xi,eta,zeta) [zeros(4,1); [xi;0]];
elseif(opt == 2)
    
    M = @(xi,eta,zeta) [zeros(4,8); [xi,zeta,0 0 xi*eta, xi*zeta,0,0; 0 0 eta, zeta, 0 0 xi*eta, eta*zeta]];
elseif(opt == 3) 
    M = @(xi,eta,zeta) [zeros(4,12); [xi 0 xi*eta 0 (1/5 - zeta^2) 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0; 0 eta 0 xi*eta 0 (1/5 - zeta^2) 0 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2) 0 xi*eta*(1/5-zeta^2)]];  
elseif(opt == 4)
    M = @(xi,eta,zeta) [zeros(2,15); ...
                       [zeros(1,12), xi, zeta zeta^3 ];...
                        zeros(1,15); ...
                       [xi 0 xi*eta 0 (1/5 - zeta^2) 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0, zeros(1,3);...
                        0 eta 0 xi*eta 0 (1/5 - zeta^2) 0 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2) 0 xi*eta*(1/5-zeta^2), zeros(1,3)]];  
elseif(opt == 5)
%     M = @(xi,eta,zeta) [[0 0 0 0 0 0 0 0 0 0 0]; [0 0 0 0 0 0 0 0 0 0 0]; [0 0 0 0 0 0 0 0 xi eta xi*eta]; [0 0 0 0 0 0 0 0 0 0 0]; [xi,zeta,0 0 xi*eta, xi*zeta,0,0 0 0 0; 0 0 eta, zeta, 0 0 xi*eta, eta*zeta 0 0 0]];
    M = @(xi,eta,zeta) [[0 0 0 0 0 0 0 0 0 0 0 0 0 0]; [0 0 0 0 0 0 0 0 0 0 0 0 0 0]; [0 0 0 0 0 0 0 0 xi eta xi*eta zeta zeta*xi zeta*eta]; [0 0 0 0 0 0 0 0 0 0 0 0 0 0]; [xi,zeta,0 0 xi*eta, xi*zeta,0,0 0 0 0 0 0 0; 0 0 eta, zeta, 0 0 xi*eta, eta*zeta 0 0 0 0 0 0]];
%     M = @(xi,eta,zeta) [eta 0 0 0 0 0]';
elseif(opt == 6)
    M = @(xi,eta,zeta) [0 0 0 0; 0 0 0 0; zeta, zeta*xi, zeta*eta, zeta*xi*eta; 0 0 0 0; 0 0 0 0;0 0 0 0];
elseif(opt == 7)
%     M = @(xi,eta,zeta) [0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 zeta, zeta*xi, zeta*eta, zeta*xi*eta; 0 0 0 0 0 0 0 0; [xi,0 xi*eta, 0 0 0 0 0; 0,eta,0,xi*eta  0 0 0 0]];
    M = @(xi,eta,zeta) [0 0 0 0 0 0 0 0 xi xi*eta 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 eta eta*xi 0 0 0; 0 0 0 0 zeta, zeta*xi, zeta*eta, zeta*xi*eta 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 xi eta xi*eta; [xi,0 xi*eta, 0 0 0 0 0  0 0 0 0 0 0 0; 0,eta,0,xi*eta  0 0 0 0 0 0 0 0 0 0 0]];
% M = @(xi,eta,zeta) [0 0 0 0 0 0 0 0 xi xi*eta 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 eta eta*xi 0 0 0 0 0 0 0; 0 0 0 0 zeta, zeta*xi, zeta*eta, zeta*xi*eta 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 xi eta xi*eta 0 0 0 0; [xi,0 xi*eta, 0 0 0 0 0  0 0 0 0 0 0 0 zeta xi*zeta 0 0; 0,eta,0,xi*eta  0 0 0 0 0 0 0 0 0 0 0 0 0 zeta eta*zeta]];
end

end