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

mkdir('Tol');
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
    app_iln_tol_all_color(dir_name, base_name, iln_min, zf_crlb_max, N_crlb_min, N_crlb_max, xf_crlb_max,  yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off, zf_max);
%     app_gpu_tol_all_color(dir_name, base_name, zf_min, iln_min, sig_min, sig_max, sig_crlb_min, sig_crlb_max, fr_sig, N_crlb_min, N_crlb_max, xf_crlb_max, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off, zf_max)
    waitbar(i/numel(finfo),w,[num2str(i), ' out of ', num2str(numel(finfo)), ' completed!']);
end
close(w);
cd(curr_path);