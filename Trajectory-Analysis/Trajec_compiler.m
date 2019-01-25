% Trajec_compiler
% This script will compile all the data from trajectory files and output a
% distribution of steps / frame for subsequent analysis
clearvars;
close all;
clc;
exp_tm = 0.04;
files = dir('*traj.mat');
dx =[];
dy =[];
dz =[];
for p = 1:numel(files)
    load(files(p).name);
    for i = 1:numel(trajec)
        ind = trajec(i).t;
        dx = [dx;q*(ncoords(ind,1)-mean(ncoords(ind,1)))];
        dy = [dy;q*(ncoords(ind,2)-mean(ncoords(ind,2)))];
        dz = [dz;q*(ncoords(ind,3)-mean(ncoords(ind,3)))];
    end
end
% histogram(D_step/(4*exp_tm));
        