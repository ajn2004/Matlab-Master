% Trajec_compiler
% This script will compile all the data from trajectory files and output a
% distribution of steps / frame for subsequent analysis
% clearvars;
close all;
clc;
exp_tm = 0.04;
% files = dir('*traj.mat');
dx =[];
dy =[];
dz =[];
% for p = 1:numel(files)
%     load(files(p).name);
    for i = 1:numel(trajec)
        ind = trajec(i).t;
        for j = 1:numel(ind)-1
            dx = [dx; q*((ncoords(ind(j),1) - ncoords(ind(j+1),1))^2)^0.5];
            dy = [dy; q*((ncoords(ind(j),2) - ncoords(ind(j+1),2))^2)^0.5];
            dz = [dz; q*((ncoords(ind(j),3) - ncoords(ind(j+1),3))^2)^0.5];
        end
        dx = [dx;q*(ncoords(ind,1)-mean(ncoords(ind,1)))];
        dy = [dy;q*(ncoords(ind,2)-mean(ncoords(ind,2)))];
        dz = [dz;q*(ncoords(ind,3)-mean(ncoords(ind,3)))];
%         dx = 
    end
% end
% histogram(D_step/(4*exp_tm));
        