%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear regression routine to identify psfs
% By using linear regression we may be able to determine with some accuracy
% wether or not an image is a PSF
%
%
%
% AJN 9/25/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% plan to get data from a mat file, or many mat files.

%% Note to developer
% figure out how to load the data. Probably be easiest to write a
% script to arrange the data before loading it into this script. The
% rest of this script will assume a data set i1 in the form
% [num_examples, mxn] = size(images) and num_examples = size(indicator)

[fname, fpath]= uigetfile('*.mat');
cd(fpath);
f_info = dir('*set.mat');
for i = 1:numel(f_info)
    load(f_info(i).name);
    if i == 1
        X = [ ones(numel(x(:,1)),1), x]; %prep data set for analysis
        Y = y;
    else
        X = vertcat(X,[ ones(numel(x(:,1)),1), x]);
        Y =vertcat(Y,y);
    end
end
int_thetas = zeros(numel(x(1,:))+1,1);  %initialize parameter set as 0s
lambda = 0.1;   %regulation factor, play around to find something useful based on trials not apart of the example set


[J , grad ] = costfunc(int_thetas, X, Y, lambda); %cost function to show first and final cost

% conv = 0.001; %convergence factor, ask J to be within 100*conv % of itself between iterations before breaking analysis
% [fin_thetas, i, Jout] = newtraphs(int_thetas, X, y.', lambda, conv); %apply netwon raphson method to determine proper values for theta parameters

options = optimset('GradObj', 'on', 'MaxIter', 1000);
    [theta] = fmincg(@(t)(costfunc(t, X, Y, lambda)), int_thetas, options);



mean(sigmoid(X*theta) > 0.5 == Y) %