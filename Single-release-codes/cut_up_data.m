%%% Cutting up the data

clearvars;
close all;
clc;
load('back_subtract.mat');
i1 = (readtiff()-mi1)/33.33;
stims = 64;
presc = 100;
imsafter = 10;
[m,n,o] = size(i1);

nstim = floor((o-presc)/stims);

%% Image cleaning
ip1 = rollingball(i1);
dip1 = [];
fms = [];
for i = 1:nstim
    ind = 100 + (i)*stims;
    dip1 = cat(3,dip1, ip1(:,:,ind:ind+imsafter-1) - ip1(:,:,ind-1));
    fms =[fms,ind:ind + imsafter-1]; 
end
