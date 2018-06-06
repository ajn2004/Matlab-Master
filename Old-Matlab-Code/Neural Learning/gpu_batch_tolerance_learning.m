%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the Batch
%
% A script that will batch apply GPU tolerances
% Just run this program
% 
% AJN 7-10-15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% User Specified Tolerances
% choose the min/max values based on the histograms of these variables
% i.e. N, off_all, xf_crlb, yf_crlb, N_crlb, off_crlb, lp, and llv
N_tol_min = 120;         % minimum number of photons
N_tol_max = 8000;       % maximum number of photons
off_min = 0;            % minimum number of offset photons
off_max = 150;           % maximum number of offset photons
N_crlb_min = 320;         % minimum error in number of photons
N_crlb_max = 8000;       % maximum error in number of photons
xf_crlb_min = 0;        % minimum error in number of x-position
xf_crlb_max = 0.05;      % maximum error in number of x-position
yf_crlb_min = 0;        % minimum error in number of y-position
yf_crlb_max = 0.05;      % maximum error in number of y-position
off_crlb_min = 0.2;     % minimum error in number of offset photons
off_crlb_max = 0.6;     % maximum error in number of offset photons
fr_unc = 0.15;
fr_unc_off = 0.08;



%% File Info (don't worry about this part)
[fname, dir_name] = uigetfile('*.mat');
curr_path = pwd;
cd(dir_name);
addpath(curr_path);
finfo = dir([dir_name,'*.mat']);
w = waitbar(0, ['0 out of ', num2str(numel(finfo)),' completed!']);
for i = 1:numel(finfo)
    base_name = finfo(i).name;
    app_gpu_tol_all_color(dir_name, base_name, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc, fr_unc_off)
    waitbar(i/numel(finfo),w,[num2str(i), ' out of ', num2str(numel(finfo)), ' completed!']);
end
close(w);
cd(curr_path);