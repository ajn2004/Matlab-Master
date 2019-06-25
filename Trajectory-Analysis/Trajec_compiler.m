% Trajec_compiler
% This script will compile all the data from trajectory files and output a
% distribution of steps / frame for subsequent analysis
% clearvars;
close all;
clc;
exp_tm = 0.04;
files = dir('*traj.mat');
dx =[];
dy =[];
dz =[];
msnr = 60;


for p = 1:numel(files)
    load(files(p).name);
    zf = func_shift_correct(ncoords(:,3)*q,framenumber,1)/q;
    ncoords(:,3) = zf(:);
    snr = fits(:,3)./(fits(:,3) + (pixw*2+1)^2*fits(:,6)).^0.5;
    for i = 1:numel(trajec)
        ind = trajec(i).t;
        indy = snr(ind) < msnr;
        ind(indy) = [];
%         for j = 1:numel(ind)-1
%             dx = [dx; q*((ncoords(ind(j),1) - ncoords(ind(j+1),1))^2)^0.5];
%             dy = [dy; q*((ncoords(ind(j),2) - ncoords(ind(j+1),2))^2)^0.5];
%             dz = [dz; q*((ncoords(ind(j),3) - ncoords(ind(j+1),3))^2)^0.5];
%         end
        if numel(ind) > 1
        dx = [dx;q*(ncoords(ind,1)-mean(ncoords(ind,1)))];
        dy = [dy;q*(ncoords(ind,2)-mean(ncoords(ind,2)))];
        dz = [dz;q*(ncoords(ind,3)-mean(ncoords(ind,3)))];
        end
%         dx = 
    end
end
% histogram(D_step/(4*exp_tm));
a = fit_hist_gauss(dz);
a(3)