clear variables

% nnn = [5 151 151; 5 171 171; 6 171 171];
nnn = [1 10 10; 1 20 20; 1 40 40; 1 60 60; 1 80 80; 1 100 100; 1 120 120];% 2 100 100];
% nnn = [5, 100 100];
for i= 1:size(nnn,1)
    
    fprintf('current IIIII: %i',i);
    disp3(i) = halfCylinderMain(nnn(i,:));
%     save('disp3','disp3')
    
end

figure
plot(disp3)

%save('convCheck_s','disp3','nnn')