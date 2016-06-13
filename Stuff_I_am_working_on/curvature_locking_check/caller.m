clear variables

% nnn = [5 151 151; 5 171 171; 6 171 171];
% nnn = [1 10 1; 1 20 2; 1 40 4; 1 60 6; 1 80 8; 1 100 10; 2 100 10];% 4 150 15];
nnn = [1 10 1; 1 20 1; 1 40 1; 1 60 1];
% nnn = [2, 100 10; 2, 120 12; 3, 120 12; 4, 120 12; 5, 170 17; 5, 250 25; 5, 300 30; 6, 300 30];
for i= 1:size(nnn,1)
    
    fprintf('current IIIII: %i',i);
    disp(i) = quarterCylinderMain(nnn(i,:));
%     save('disp3','disp3')
    
end

totnel = prod(nnn,2);

% figure
plot(prod(nnn,2),disp)
ylabel('Defl'); xlabel('Number of elements')
% save('zSE4A3','totnel','disp')



