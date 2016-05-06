function [M] = getInterPolMatrix(opt)

if(opt == 1)
    M = @(xi,eta,zeta) [[0 0 0 0]; [0 0 0 0]; [0 0 0 0] ;zeros(1,4); [xi,0 xi*eta, 0; 0,eta,0,xi*eta]];
elseif(opt == 2)  
    M = @(xi,eta,zeta) [zeros(4,8); [xi,zeta,0 0 xi*eta, xi*zeta,0,0; 0 0 eta, zeta, 0 0 xi*eta, eta*zeta]];
elseif(opt == 3) 
    M = @(xi,eta,zeta) [zeros(4,12); [xi 0 xi*eta 0 (1/5 - zeta^2) 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0; 0 eta 0 xi*eta 0 (1/5 - zeta^2) 0 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2) 0 xi*eta*(1/5-zeta^2)]];  
elseif(opt == 4)
    %GOODWITHANS--
    M =  @(x,y,z)reshape([x,0.0,0.0,0.0,0.0,0.0,x.*y,0.0,0.0,0.0,0.0,0.0,0.0,y,0.0,0.0,0.0,0.0,0.0,x.*y,0.0,0.0,0.0,0.0,0.0,0.0,z,0.0,0.0,0.0,0.0,0.0,x.*z,0.0,0.0,0.0,0.0,0.0,y.*z,0.0,0.0,0.0,0.0,0.0,x.*y.*z,0.0,0.0,0.0,0.0,0.0,0.0,x,0.0,0.0,0.0,0.0,0.0,y,0.0,0.0,0.0,0.0,0.0,x.*y,0.0,0.0],[6,11]);
    %GOODWITHANS--
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