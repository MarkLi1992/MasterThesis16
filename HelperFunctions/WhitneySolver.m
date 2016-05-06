function WhitneySolver()
a = 1; b = 1; zz = 0.001;
p = -1000;
%  
%         EL = 50E9;    
%         ET = 9E9;   
%         nuLT = 0.22;    
%         GLT = 5E9; 
%         GTT = 3.2E9;     
%         nuTL = ET/EL*nuLT;
      
        EL = 174.6E9;  ET = 7E9;  nuLT = 0.25;    GLT = 3.5E9;   GTT = 1.4E9; nuTL = ET/EL*nuLT;
        
QLT = [EL/(1-nuLT*nuTL) (nuLT*ET)/(1-nuLT*nuTL) 0;...
       (nuLT*ET)/(1-nuLT*nuTL) ET/(1-nuLT*nuTL) 0;...
       0 0 GLT];

QLTtilde = [GLT 0; 0 GTT];

ang = [15 -15]*pi/180;
coord = linspace(-zz/2, zz/2,length(ang)+1); 

[maxabs_a] = WHITNEY(a,b,p,QLT,QLTtilde, ang, coord)
fprintf('Max defl: %d\n' ,maxabs_a)
end

function [maxabs_a] = WHITNEY(a,b,q0,QLT,QLTtilde, ang, coord)
% Example Input can be found furthest down in the code

[A,B,D,Atilde] = ABDAtilde(QLT,QLTtilde, ang, coord);

R = a/b;

A11 = A(1,1); A12 = A(1,2); A66 = A(3,3); A22 = A(2,2);
D11 = D(1,1); D22 = D(2,2); D12 = D(1,2); D66 = D(3,3);
B16 = B(1,3); B26 = B(2,3); 

xx = a/2; yy = b/2;
w = 0;
toM = 100; toN = 100;
for m = 1:2:toM
    for n = 1:2:toN
        
        Hmn = (((A11*m^2 + A66*n^2*R^2)*(A66*m^2 + A22*n^2*R^2) - ...
               (A12 + A66)^2*m^2*n^2*R^2)*(D11*m^4 + 2*(D12 + 2*D66)*m^2*n^2*R^2 + ...
               D22*n^4*R^4) + 2*m^2*n^2*R^2*(A12 + A66)*(3*B16*m^2 + ...
               B26*n^2*R^2)*(B16*m^2 + 3*B26*n^2*R^2) - ...
               n^2*R^2*(A66*m^2 + A22*n^2*R^2)*(3*B16*m^2 + B26*n^2*R^2)^2 - ...
               m^2*(A11*m^2 + A66*n^2*R^2)*(B16*m^2 + 3*B26*n^2*R^2)^2);

        qmn = 16*q0/pi^2/m/n; if(mod(m,2) == 0 || mod(n,2) == 0); qmn = 0; end;     

        Gmn = (qmn*R^4*b^4)/(pi^4*Hmn)*((A11*m^2 + A66*n^2*R^2)*(A66*m^2 + A22*n^2*R^2) - (A12 + A66)^2*m^2*n^2*R^2);   
        
        w = w + Gmn*sin(m*pi*xx/a)*sin(n*pi*yy/b);
        
    end
end

maxabs_a = w;

end
function [A,B,D,Atilde] = ABDAtilde(QLT,QLTtilde, ang, coord)
% Function calculating A, B and D matrix for a laminate

% Inputs:
%        QLT: Local stiffness matrix
%        ang: vector with the angles of the lamina (starting from min. z)
%        coord: vector with z-coordinates of the lamina interfaces
%        alfaLT: Local heat expansion vector
%        deltaT: Difference in temperature

%    Elias Börjesson, Fredrik Ekre

T1 = @(theta) [cos(theta)^2 sin(theta)^2 2*sin(theta)*cos(theta);...
    sin(theta)^2 cos(theta)^2 -2*sin(theta)*cos(theta);...
    -sin(theta)*cos(theta) sin(theta)*cos(theta) (cos(theta)^2-sin(theta)^2)];

T2 = @(theta) [cos(theta)^2 sin(theta)^2 sin(theta)*cos(theta);...
    sin(theta)^2 cos(theta)^2 -sin(theta)*cos(theta);...
    -2*sin(theta)*cos(theta) 2*sin(theta)*cos(theta) (cos(theta)^2-sin(theta)^2)];

T1tilde = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

A = zeros(3); B = A; D = A; Atilde = zeros(2);
Qxy = @(theta) T1(theta)^-1*QLT*T2(theta);
Qxytilde = @(theta) T1tilde(theta)^-1*QLTtilde*T1tilde(theta);
nla = length(ang);

for ii = 1:nla
    A = A + 1/1*Qxy(ang(ii))*(coord(ii+1)^1-coord(ii)^1);
    B = B + 1/2*Qxy(ang(ii))*(coord(ii+1)^2-coord(ii)^2);
    D = D + 1/3*Qxy(ang(ii))*(coord(ii+1)^3-coord(ii)^3);
    Atilde = Atilde + Qxytilde(ang(ii))*(coord(ii+1)-coord(ii));
end
end

