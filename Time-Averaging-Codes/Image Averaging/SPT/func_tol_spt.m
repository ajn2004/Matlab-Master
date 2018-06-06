function func_tol_spt(fname,q, llv_max, llv_min, sigma2_min, sigma2_max, sig_crlb_min, sig_crlb_max, fr_unc_sig, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU Tolerances
% 
% This script will apply tolerances to the GPU localized data and does so
% based off of the CRLB values obtained for a given data set
%
% Written by AJN 7-10-15 Based off Hess Lab's Apply Tolerance program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fname);  % load file
%% assign variables to a temporary variable
N_temp = N_all;
N_crlb_temp = N_crlb;
xf_all_temp = xf_all;
xf_crlb_temp = xf_crlb;
yf_all_temp = yf_all;
yf_crlb_temp = yf_crlb;
off_all_temp = -off_all;
off_crlb_temp = -off_crlb;
llv_temp = llv;
sigx_temp = sigx_all;
sigy_temp = sigy_all;
sigx_crlb_temp = sigx_crlb;
sigy_crlb_temp = sigy_crlb;
framenum_all_temp = framenum_all;
fr_N = real(N_crlb.^0.5./N_all);
fr_off = real(off_crlb.^0.5./off_all);
fr_sigx = real(sigx_crlb_temp.^0.5./sigx_temp);
fr_sigy = real(sigy_crlb_temp.^0.5./sigy_temp);
sig_max = sigma2_max/(2*q*1000);
sig_min = sigma2_min/(2*q*1000);
num_temp = number;

%% load multicolor info if it exists
if exist('nrat')
    nrat_temp = nrat;
    ratio_temp = ratio;
    clear nrat ratio
end
%% clear variables
clear N_all N_crlb xf_all xf_crlb yf_all yf_crlb off_all off_crlb llv lp framenum_all sigx_all sigy_all sigx_crlb sigy_crlb number
%% determine index of acceptable molecules
% llv_temp <= llv_max & llv_temp >= llv_min & ... llv temps
index = find(sigx_crlb_temp >= sig_crlb_min & sigy_crlb_temp >= sig_crlb_min & sigx_crlb_temp <= sig_crlb_max & sigy_crlb_temp <= sig_crlb_max &... sigma _crlb tols
    sigx_temp <= sig_max & sigy_temp <= sig_max & sigx_temp >= sig_min & sigy_temp >= sig_min & ... sigma tols
    fr_sigx <= fr_unc_sig & fr_sigy <= fr_unc_sig & fr_N <= fr_unc_N & fr_off <= fr_unc_off & ...  frac unc tols
    N_temp <= N_tol_max & N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max & ... % N tols
    xf_crlb_temp <= xf_crlb_max & xf_crlb_temp >= xf_crlb_min & yf_crlb_temp <= yf_crlb_max & yf_crlb_temp >= yf_crlb_min & ... coord tols
    off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max & off_all_temp >= off_min & off_all_temp <= off_max); % offset tols

%% redefine variables to be only of acceptable molecules
N = N_temp(index);
N_crlb = N_crlb_temp(index);
xf_all = xf_all_temp(index);
xf_crlb = xf_crlb_temp(index);
yf_all = yf_all_temp(index);
yf_crlb = yf_crlb_temp(index);
off_all = off_all_temp(index);
off_crlb = off_crlb_temp(index);
% llv = llv_temp(index);
% lp = lp_temp(index);
sigx_all = sigx_temp(index);
sigy_all = sigy_temp(index);
sigx_crlb = sigx_crlb_temp(index);
sigy_crlb = sigy_crlb_temp(index);
framenum_all = framenum_all_temp(index);
number = num_temp(index);

if exist('nrat_temp')  % multicolor if it exists
    nrat = nrat_temp(index);
    ratio = ratio_temp(index);
    clear nrat_temp ratio_temp
end
total_molecules = numel(N);
%% clear temporary variables
clear num_temp N_temp N_crlb_temp xf_all_temp xf_crlb_temp yf_all_temp yf_crlb_temp off_all_temp off_crlb_temp llv_temp lp_temp framenum_all_temp sigx_temp sigy_temp sigy_crlb_temp sigx_crlb_temp

%% save workspace and append '_tol' to file name
save([fname(1:end-4), '_tol']);