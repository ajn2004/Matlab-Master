clearvars; close all; clc;
% Z-Unc Sim
% This simulation takes in values for sigx, sigy, and the uncertainty in
% those measurements to return a simulated distribution in Z
simsize = 100000;
q = 0.128;
sx = 1.4564;
sy = 2.1437;
su = 0.1280;

dsx = normrnd(sx,su,simsize,1);
dsy = normrnd(sy,su,simsize,1);

load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Single-release-codes\z_calib.mat');
zfs = getdz(dsx,dsy,cal.z_cal);
histogram(zfs)
xlabel('Estimated Z Position nm')
ylabel('Frequency')
std(zfs)