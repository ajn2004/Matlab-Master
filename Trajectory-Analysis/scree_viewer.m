% Quick Look
%
% This will quickly organize the trajectory data to show number of
% trajectories as a function of dmax

clearvars; close all; clc;
files = dir('*traj.mat');

trajs = zeros(numel(files),1);
trajl = [];
dms = trajs;
for i = 1:numel(files)
    load(files(i).name,'trajec','dmax','total_molecules');
    trajs(i) = numel(trajec);
    trajl(i) = 0;
    for j = 1:trajs(i)
        trajl(i) = trajl(i) + numel(trajec(j).t);
    end
    trajl(i) = trajl(i)/trajs(i);
    dms(i) = dmax;
    trajp(i) = trajs(i)./total_molecules;
    clear trajec dmax
end

plot(dms,trajs,'.');
figure
plot(dms,trajl,'.');