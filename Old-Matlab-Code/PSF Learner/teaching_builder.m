%  This script will build data necessary to teach a computer how to
%  identify point spread functions
%
% ajn 9/28/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clear all
% clc
% 
% [fname, fpath] = uigetfile('*tol.mat','Grab a file that has the grab images');
% addpath(pwd);
% cd(fpath);
% file_inf = dir('*.mat');

for l = 1:numel(file_inf)
    l
    num_samp = 1000;
    load(file_inf(l).name, 'i5');
    m = randperm(numel(i5(1,1,:)));
    if 2*num_samp > numel(i5(1,1,:));
        num_samp = floor(numel(i5(1,1,:))/2);
    end
    for i = 1:2*num_samp; % create 2 times as much positive data to make equal with number of negative data
        fprintf('%d%% done\n',round(100*(i)/(num_samp*3)));
        ease = i5(:,:,m(i));
        x(i,:) = ease(:).';
        y(i,1) = 0;
    end
%     n = randperm(numel(grab_failed(1,1,:)));
%     for j = 1:num_samp;
%         fprintf('%d%% done\n',round(100*(i+j)/(num_samp*3)));
%         ease = i5(:,:,n(j));
%         x(i+j,:) = ease(:).';
%         y(i+j,1) = 0;
%     end
%     mean_val = mean(x(:));
%     i1 = ones(numel(x(1,:))^0.5)*mean_val; 
%     for k = 1: num_samp
%     i2 = double(imnoise(uint16(i1), 'poisson'));
%     x = vertcat(x,i2(:).');
%     y = vertcat(y,0);
%     end
    save([file_inf(l).name(1:end-4), '_set.mat'], 'x','y')
    clear x y grab_failed grab_passed
end