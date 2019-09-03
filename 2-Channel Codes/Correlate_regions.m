%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 Color transform
% 
% My attempt at making a 2-Image Region correlator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all
clc;

fname = 'scan_at_10nm_2.tif';
fpath = 'C:\Users\AJN Lab\Dropbox\Data\7-15-19 2-Color Configs\double_scan\';
i1 = readtiff([fpath,fname]);
[m,n,o] = size(i1);
cutn = round(n/2);
im1 = i1(:,1:cutn,:);
im2 = i1(:,cutn+1:end,:);
imagesc(mean(im1,3)+mean(im2,3));