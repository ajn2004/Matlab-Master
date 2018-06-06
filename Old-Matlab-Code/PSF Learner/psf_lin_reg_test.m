%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PSF lin reg tester
%
%
% a script to test the results of the PSF linear regulator
%
% AJN 9/29/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

load('learned_theta.mat');
[fname, fpath] = uigetfile('*set.mat','Grab a test file');

cd(fpath);
file_inf = dir('*set.mat');
for k = 1:101
    k
for i = 1:numel(file_inf)
    load(file_inf(i).name);
    
    X = [ones(numel(x(:,1)),1), double(x)];
    
    ave(i,k) = mean((sigmoid(X*theta)  > (k-1)/100) == y);
    yes(i,k) = mean((sigmoid(X(y==1, :)*theta)  > (k-1)/100) == y(y==1));
    nos(i,k) = mean((sigmoid(X(y==0, :)*theta)  > (k-1)/100) == y(y==0));
end
end

plot(mean(ave));
hold on
plot(mean(yes),'g');
plot(mean(nos),'r');
legend('Probability','Yes','No')
