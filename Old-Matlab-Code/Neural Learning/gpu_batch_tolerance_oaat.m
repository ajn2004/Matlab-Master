%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the Batch
%
% A script that will batch apply GPU tolerances
% Just run this program This program will only analyze the selected file
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
N_crlb_min = 200;         % minimum error in number of photons
N_crlb_max = 1800;       % maximum error in number of photons
xf_crlb_min = 0;        % minimum error in number of x-position
xf_crlb_max = 0.03;      % maximum error in number of x-position
yf_crlb_min = 0;        % minimum error in number of y-position
yf_crlb_max = 0.03;      % maximum error in number of y-position
off_crlb_min = 0.13;     % minimum error in number of offset photons
off_crlb_max = 1.35;     % maximum error in number of offset photons
N_tol_min = 90;         % minimum number of photons
N_tol_max = 1100;       % maximum number of photons
% llv_min = -3*10^15;     % minimum log likelihood value
% llv_max = 0;            % maximum log likelihood value
off_min = 4.8;            % minimum number of offset photons
off_max = 50;           % maximum number of offset photons
lp_min = 0;             % minimum error in loc uncertainty
lp_max = 1;             % maximum error in loc uncertainty



%% File Info (don't worry about this part)
[fname, dir_name] = uigetfile('*.mat');
curr_path = pwd;
cd(dir_name);
addpath(curr_path);
base_name = fname;
app_gpu_tol_all_color(dir_name, base_name, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, lp_min, lp_max, off_min, off_max)
cd(curr_path);