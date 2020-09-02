
% calling kernfun2 which should generate sigmaList kernels quicker

clear all; close all; clc;

data = loaddata();

sigmaList = -8:0.2:6;
distmat0 = kerfun2(data, data, 10.^sigmaList, 0);
distmat1 = kerfun2(data, data, 10.^sigmaList, 1);
for i = 1:length(sigmaList)
    s = sigmaList(i);
    csvwrite(strcat('reviewKernels/dim_0_1e',num2str(round(s/.1)*.1),'.csv'),distmat0(:,:,i));
    csvwrite(strcat('reviewKernels/dim_1_1e',num2str(round(s/.1)*.1),'.csv'),distmat1(:,:,i));
end

