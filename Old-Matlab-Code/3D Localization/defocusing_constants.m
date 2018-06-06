%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blank Document
%
% This is the template that I always use for making any kind of script
%
%
% AJN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

focal_plane = 550; % frame number focal plane is found on
[fname, fpath] = uigetfile('*.mat');
cd(fpath);
files = dir('*dast.mat');

sigxa = [];
sigya = [];
fnuma = [];
for i = 1:numel(files)
    load(files(i).name);
%     ind = y == 1;
%     sum(y)
%     numel(sigx_all)
    fnuma = [fnuma;framenum_all];
    sigxa = [sigxa;sigx_all];
%     sigxc = [sigxc;sigx_crlb];
    sigya = [sigya;sigy_all];
%     sigyc = [sigyc;sigy_crlb];
%     iloca = [iloca,iloc(:,ind)];
    clearvars -except files sigxa sigya iloca fnuma i sigxc sigyc focal_plane
end
z = (fnuma - focal_plane)*2;

