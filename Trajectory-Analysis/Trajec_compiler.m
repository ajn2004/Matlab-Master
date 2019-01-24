% Trajec_compiler
% This script will compile all the data from trajectory files and output a
% distribution of steps / frame for subsequent analysis
clearvars;
close all;
clc;
exp_tm = 0.04;
files = dir('*traj.mat');
D_step =[];
for p = 1:numel(files)
    load(files(p).name);
    for i = 1:numel(trajec)
        for k = 1:numel(trajec(i).t)-1
            id = trajec(i).t(k);
            id2 = trajec(i).t(k+1);
            D_step(numel(D_step)+1) = q^2*((ncoords(id,1) - ncoords(id2,1))^2 +(ncoords(id,2) - ncoords(id2,2))^2 + (ncoords(id,3) - ncoords(id2,3))^2);
        end
    end
end
histogram(D_step/(4*exp_tm));
        