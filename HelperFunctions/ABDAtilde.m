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