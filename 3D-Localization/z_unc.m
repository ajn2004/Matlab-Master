function [unc] = z_unc(sx, sy, ux, uy)

simsize = 100;
dsx = normrnd(sx,ux,simsize,1);
dsy = normrnd(sy,uy,simsize,1);
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Single-release-codes\z_calib.mat');
zfs = getdz(dsx,dsy,cal.z_cal);
unc = std(zfs);