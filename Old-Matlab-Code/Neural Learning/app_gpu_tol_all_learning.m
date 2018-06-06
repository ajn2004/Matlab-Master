function app_gpu_tol_all_color_learning(dir_name, base_name, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc, fr_unc_off)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU Tolerances
% 
% This script will apply tolerances to the GPU localized data and does so
% based off of the CRLB values obtained for a given data set
%
% Written by AJN 7-10-15 Based off Hess Lab's Apply Tolerance program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dir_name, base_name]);  % load file
%% assign variables to a temporary variable
N_temp = N;
N_crlb_temp = N_crlb;
xf_all_temp = xf_all;
xf_crlb_temp = xf_crlb;
yf_all_temp = yf_all;
yf_crlb_temp = yf_crlb;
off_all_temp = off_all;
off_crlb_temp = off_crlb;
llv_temp = llv;
lp_temp  = lp;
framenum_all_temp = framenum_all;
fr_N = N_crlb.^0.5./N;
fr_off = off_crlb.^0.5./off_all;

%% load multicolor info if it exists
if exist('nrat');
    nrat_temp = nrat;
    ratio_temp = ratio;
    clear nrat ratio
end
%% clear variables
clear N N_crlb xf_all xf_crlb yf_all yf_crlb off_all off_crlb llv lp framenum_all
%% determine index of acceptable molecules
index = find (fr_N < fr_unc & fr_off < fr_unc_off & N_temp < N_tol_max & N_temp > N_tol_min & N_crlb_temp > N_crlb_min & N_crlb_temp < N_crlb_max & xf_crlb_temp < xf_crlb_max & xf_crlb_temp > xf_crlb_min & yf_crlb_temp < yf_crlb_max & yf_crlb_temp > yf_crlb_min & off_crlb_temp > off_crlb_min & off_crlb_temp < off_crlb_max & off_all_temp >off_min & off_all_temp <off_max);

%% redefine variables to be only of acceptable molecules
N = N_temp(index);
N_crlb = N_crlb_temp(index);
xf_all = xf_all_temp(index);
xf_crlb = xf_crlb_temp(index);
yf_all = yf_all_temp(index);
yf_crlb = yf_crlb_temp(index);
off_all = off_all_temp(index);
off_crlb = off_crlb_temp(index);
llv = llv_temp(index);
lp = lp_temp(index);
framenum_all = framenum_all_temp(index);
if exist('nrat_temp');  % multicolor if it exists
    nrat = nrat_temp(index);
    ratio = ratio_temp(index);
    clear nrat_temp ratio_temp
end
total_molecules = numel(N);
%% clear temporary variables
clear N_temp N_crlb_temp xf_all_temp xf_crlb_temp yf_all_temp yf_crlb_temp off_all_temp off_crlb_temp llv_temp lp_temp framenum_all_temp

%% save workspace and append '_tol' to file name
save([base_name(1:end-4), '_tol']);