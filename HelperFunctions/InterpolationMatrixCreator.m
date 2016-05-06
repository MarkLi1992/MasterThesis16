clear all;

syms x y z
s(1).terms = {};
s(2).terms = {};
s(3).terms = {z, z*x , z*y, x*y*z};
s(4).terms = {};
s(5).terms = {x, z, x*y, x*z}%{x, x*y, (1/5 - z^2), x*(1/5 - z^2), y*(1/5 - z^2), x*y*(1/5 - z^2)};
s(6).terms = {y, z, x*y, y*z}%{y, x*y, (1/5 - z^2), x*(1/5 - z^2), y*(1/5 - z^2), x*y*(1/5 - z^2)};

n = size([s.terms],2);

Mtmp = sym(zeros(6,n));
itr = 1;
for i=1:6
   
    for j = 1:size(s(i).terms,2)
       
        Mtmp(i,itr) = s(i).terms(j);
        itr = itr + 1;
    end
    
end
Mtmp
matlabFunction(Mtmp)

% [xi 0 xi*eta 0 (1/5 - zeta^2) 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0; 0 eta 0 xi*eta 0 (1/5 - zeta^2) 0 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2) 0 xi*eta*(1/5-zeta^2)]]