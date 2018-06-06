%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch trajectories
%
% a trajectory code to batch run on a data set
%
% ajn 2/6/16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

% dmax = 20;

i = 3;
files = dir('*.mat');
for dmax = 5:5:500
    for i = 1:numel(files)
        func_traj_analysis(files(i).name,dmax);
       
    end
end

