
% a1 = [0 90 0]; a2 = [0 90 0 90 90 0 90 0]; a3 = [0 90 0 90 0 90 0 90 90 0 90 0 90 0 90 0];
a1 = [15 -15]; a2 = [30 -30]; a3 = [45 -45];
setupinfo(1).funk = @Solid8AnsEasSR;
setupinfo(1).M = 4;
setupinfo(1).angles = a1;

setupinfo(2).funk = @Solid8AnsEasSR;
setupinfo(2).M = 4;
setupinfo(2).angles = a2;

setupinfo(3).funk = @Solid8AnsEasSR;
setupinfo(3).M = 4;
setupinfo(3).angles = a3;
%---
setupinfo(4).funk = @Solid8StressRecLayered;
setupinfo(4).M = 1;
setupinfo(4).angles = a1;

setupinfo(5).funk = @Solid8StressRecLayered;
setupinfo(5).M = 1;
setupinfo(5).angles = a2;

setupinfo(6).funk = @Solid8StressRecLayered;
setupinfo(6).M = 1;
setupinfo(6).angles = a3;
%----
setupinfo(7).funk = @Solid8StressRecLayered;
setupinfo(7).M = 2;
setupinfo(7).angles = a1;

setupinfo(8).funk = @Solid8StressRecLayered;
setupinfo(8).M = 2;
setupinfo(8).angles = a2;

setupinfo(9).funk = @Solid8StressRecLayered;
setupinfo(9).M = 2;
setupinfo(9).angles = a3;



for i=1:9
    
    maxabs_a(i) = LSF_main(setupinfo(i));
    
end

