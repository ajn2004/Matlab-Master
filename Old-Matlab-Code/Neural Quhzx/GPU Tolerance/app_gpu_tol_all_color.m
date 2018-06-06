function app_gpu_tol_all_color(dir_name, base_name, zf_min,iln_min, sig_min, sig_max, sig_crlb_min, sig_crlb_max, fr_sig, N_crlb_min, N_crlb_max, xf_crlb_max, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off, zf_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU Tolerances
% 
% This script will apply tolerances to the GPU localized data and does so
% based off of the CRLB values obtained for a given data set
%
% Written by AJN 7-10-15 Based off Hess Lab's Apply Tolerance program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bname = base_name;
load([dir_name, base_name]);  % load file
%% assign variables to a temporary variable
zf_all = getdz(sigx_all,sigy_all)/q;
N_temp = N;
N_crlb_temp = N_crlb;
xf_all_temp = xf_all;
xf_crlb_temp = xf_crlb;
yf_all_temp = yf_all;
yf_crlb_temp = yf_crlb;
off_all_temp = off_all;
off_crlb_temp = off_crlb;
llv_temp = llv;
sigx_temp = sigx_all;
sigy_temp = sigy_all;
sigx_crlb_temp = sigx_crlb;
sigy_crlb_temp = sigy_crlb;
framenum_all_temp = framenum_all;
iln = llv./N;
fr_N_t = N_crlb.^0.5./N;
fr_off_t = off_crlb.^0.5./off_all;
fr_sigx_t = sigx_crlb_temp.^0.5./sigx_temp;
fr_sigy_t = sigy_crlb_temp.^0.5./sigy_temp;
% sig_max = sig_max/(2*q*1000);
% sig_min = sigma2_min/(2*q*1000);
iloc_temp = iloc(:,logical(y));
zf_nm_temp = zf_all;
%% load multicolor info if it exists
if exist('nrat')
    nrat_temp = nrat;
    ratio_temp = ratio;
    clear nrat ratio
end
%% clear variables
clear N iloc y N_crlb xf_all xf_crlb yf_all yf_crlb off_all off_crlb llv lp framenum_all sigx_all sigy_all sigx_crlb sigy_crlb zf_nm
%% determine index of acceptable molecules
% index = find(xf_all_temp >= 5 & yf_all_temp >= 5 & xf_all_temp <= max(xf_all_temp) - 5 & yf_all_temp <= max(yf_all_temp) - 5 & iln >= iln_min & zf_nm_temp <= zf_max & zf_nm_temp >= zf_min & ... llv temps
%     sigx_crlb_temp >= sig_crlb_min & sigy_crlb_temp >= sig_crlb_min & sigx_crlb_temp <= sig_crlb_max & sigy_crlb_temp <= sig_crlb_max &... sigma _crlb tols
%     sigx_temp <= sig_max & sigy_temp <= sig_max & sigx_temp >= sig_min & sigy_temp >= sig_min & ... sigma tols
%     fr_sigx_t <= fr_sig & fr_sigy_t <= fr_sig & fr_N_t <= fr_unc_N & fr_off_t <= fr_unc_off & ...  frac unc tols
%     N_temp <= N_tol_max & N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max & ... % N tols
%     xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max & ... coord tols
%     off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max & off_all_temp >= off_min & off_all_temp <= off_max); % offset tols

index = find(xf_all_temp >= 5 & yf_all_temp >= 5 &  iln >= iln_min & zf_nm_temp <= zf_max & zf_nm_temp >= zf_min & ... llv temps
    sigx_crlb_temp >= sig_crlb_min & sigy_crlb_temp >= sig_crlb_min & sigx_crlb_temp <= sig_crlb_max & sigy_crlb_temp <= sig_crlb_max &... sigma _crlb tols
    sigx_temp <= sig_max & sigy_temp <= sig_max & sigx_temp >= sig_min & sigy_temp >= sig_min & ... sigma tols
    fr_sigx_t <= fr_sig & fr_sigy_t <= fr_sig & fr_N_t <= fr_unc_N & fr_off_t <= fr_unc_off & ...  frac unc tols
    N_temp <= N_tol_max & N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max & ... % N tols
    xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max & ... coord tols
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
llv = llv_temp(index);
fr_N = fr_N_t(index);
fr_off = fr_off_t(index);
fr_sigx = fr_sigx_t(index);
fr_sigy = fr_sigy_t(index);
% lp = lp_temp(index);
sigx_all = sigx_temp(index);
sigy_all = sigy_temp(index);
sigx_crlb = sigx_crlb_temp(index);
sigy_crlb = sigy_crlb_temp(index);
framenum_all = framenum_all_temp(index);
try
iloc = iloc_temp(:,index);
catch lsterr
end
zf_all = zf_nm_temp(index);
% zf_nm = getdz(sigx_all,sigy_all);
if exist('nrat_temp')  % multicolor if it exists
    nrat = nrat_temp(index);
    ratio = ratio_temp(index);
    clear nrat_temp ratio_temp
end
total_molecules = numel(N);
%% clear temporary variables
clear iloc_temp sig_max sig_min sigma2_max sigma2_min xf_crlb_max xf_crlb_min yf_crlb_max yf_crlb_min llv_min llv_max fr_unc_N fr_unc_off fr_unc_sig N_crlb_max N_crlb_min N_tol_max N_tol_min off_crlb_max off_crlb_min off_max off_min sig_crlb_max sig_crlb_min index fr_N_t fr_off_t fr_sigx_t fr_sigy_t N_temp N_crlb_temp xf_all_temp xf_crlb_temp yf_all_temp yf_crlb_temp off_all_temp off_crlb_temp llv_temp lp_temp framenum_all_temp sigx_temp sigy_temp sigy_crlb_temp sigx_crlb_temp
[pwd,'\Tol\',bname(1:end-4), '_tol']
%% save workspace and append '_tol' to file name
save([pwd,'\Tol\',bname(1:end-4), '_tol']);