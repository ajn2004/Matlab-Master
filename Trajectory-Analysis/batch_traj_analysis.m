%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch trajectories
%
% a trajectory code to batch run on a data set
%
% ajn 2/6/16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

dmax = 60;

% i = 3;
files = dir('*.mat');
% ds = 10.^([0 0.5 1 1.5 2 2.5 3 3.5 4]);
% for dmax = 10:10:100
    for i = 1:numel(files)
        func_traj_analysis(files(i).name,dmax);  
    end
% end

