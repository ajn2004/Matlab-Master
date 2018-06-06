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

mkdir('Tol1');
% if ~isempty('Tol')
% sendit2('Tol');
% end
%% User Specified Tolerances
% choose the min/max values based on the histograms of these variables
% i.e. N, off_all, xf_crlb, yf_crlb, N_crlb, off_crlb, lp, and llv
% Tolerance Variables
 % Tolerance Variables
%    N_tol_min = 20;        % minimum number of photons
%    N_tol_max = 10000;       % maximum number of photons
%      off_min = 0;          % minimum number of offset photons
%      off_max = 5;         % maximum number of offset photons
%   N_crlb_min = 0;          % minimum variance in number of photons
%   N_crlb_max = 10000;      % maximum variance in number of photons
%   sigma2_min = 200;        % minimum sigma value
%   sigma2_max = 800;        % maximum sigma value
%  xf_crlb_min = 0;          % minimum error in number of x-position
%  xf_crlb_max = 0.16;       % maximum error in number of x-position
%  yf_crlb_min = 0;          % minimum error in number of y-position
%  yf_crlb_max = 0.16;       % maximum error in number of y-position
% off_crlb_min = 0.0;          % minimum error in number of offset photons
% off_crlb_max = 5;          % maximum error in number of offset photons
% sig_crlb_min = 0.0;
% sig_crlb_max = 0.2;
     llv_max = 0;
     llv_min = -80000;
%     fr_unc_N = 0.4;        % fractional uncertainty in N
%   fr_unc_off = 0.15;        % Fractional uncertainty in offset
%   fr_unc_sig = 0.175;       % Fractional uncertainty in width
load('these_tol.mat');
%% File Info (don't worry about this part)
[fname, dir_name] = uigetfile('*.mat');
curr_path = pwd;
cd(dir_name);
addpath(curr_path);
finfo = dir([dir_name,'*.mat']);
w = waitbar(0, ['0 out of ', num2str(numel(finfo)),' completed!']);
for i = 1:numel(finfo)
    base_name = finfo(i).name;
    app_iln_tol_all_color(dir_name, base_name, llv_max, llv_min, sigma2_min, sigma2_max, sig_crlb_min, sig_crlb_max, fr_unc_sig, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off, zf_max);
    waitbar(i/numel(finfo),w,[num2str(i), ' out of ', num2str(numel(finfo)), ' completed!']);
end
close(w);
cd(curr_path);