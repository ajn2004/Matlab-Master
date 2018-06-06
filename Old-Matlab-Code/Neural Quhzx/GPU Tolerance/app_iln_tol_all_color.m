function app_iln_tol_all_color(dir_name, base_name, iln_min, zf_crlb_max, N_crlb_min, N_crlb_max, xf_crlb_max,  yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off, zf_max);
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
iln = llv./N;
fr_N = N_crlb.^0.5./N;
fr_off = off_crlb.^0.5./off_all;
N_temp = N;
N_crlb_temp = N_crlb;
xf_all_temp = xf_all;
xf_crlb_temp = xf_crlb;
yf_all_temp = yf_all;
yf_crlb_temp = yf_crlb;
zf_all_temp = zf_all;
zf_crlb_temp = zf_crlb;
off_all_temp = off_all;
off_crlb_temp = off_crlb;
llv_temp = llv;
framenum_all_temp = framenum_all;
% iloc_temp = iloc(:,logical(y));

%% load multicolor info if it exists
if exist('nrat');
    nrat_temp = nrat;
    ratio_temp = ratio;
    clear nrat ratio
end
%% clear variables
clear N iloc y N_crlb xf_all xf_crlb yf_all yf_crlb off_all off_crlb llv lp framenum_all sigx_all sigy_all sigx_crlb sigy_crlb zf_all zf_crlb
%% determine index of acceptable molecules
% index = find(xf_all_temp >= 2 & yf_all_temp >= 2 & xf_all_temp <= max(xf_all_temp) - 2 &...
%     yf_all_temp <= max(yf_all_temp) - 2 & iln >= iln_min & abs(zf_all_temp) <= zf_max &...
%     fr_N <= fr_unc_N & fr_off <= fr_unc_off & N_temp <= N_tol_max &...
%     N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max &...
%     xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max &...
%     off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max &...
%     off_all_temp >= off_min & off_all_temp <= off_max);

% Uncomment this section when you solve the bugs with Z-fitting  xf_all_temp <= max(xf_all_temp) - 2 &yf_all_temp <= max(yf_all_temp) - 2 & 
index = find(xf_all_temp >= 2 & yf_all_temp >= 2 &...
    iln >= iln_min & abs(zf_all_temp) <= zf_max &...
    fr_N <= fr_unc_N & fr_off <= fr_unc_off & N_temp <= N_tol_max &...
    N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max &...
    xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max &...
    zf_crlb_temp <= zf_crlb_max & off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max &...
    off_all_temp >= off_min & off_all_temp <= off_max);

%% redefine variables to be only of acceptable molecules
N = N_temp(index);
N_crlb = N_crlb_temp(index);
xf_all = xf_all_temp(index);
xf_crlb = xf_crlb_temp(index);
yf_all = yf_all_temp(index);
yf_crlb = yf_crlb_temp(index);
zf_all = zf_all_temp(index);
zf_crlb = zf_crlb_temp(index);
off_all = off_all_temp(index);
off_crlb = off_crlb_temp(index);
llv = llv_temp(index);


framenum_all = framenum_all_temp(index);
try
iloc = iloc_temp(:,index);
catch lsterr
end

% zf_nm = getdz(sigx_all,sigy_all);
if exist('nrat_temp')  % multicolor if it exists
    nrat = nrat_temp(index);
    ratio = ratio_temp(index);
    clear nrat_temp ratio_temp
end
total_molecules = numel(N);
%% clear temporary variables
clear zf_crlb_temp zf_all_temp iloc_temp sig_max sig_min sigma2_max sigma2_min xf_crlb_max xf_crlb_min yf_crlb_max yf_crlb_min llv_min llv_max fr_unc_N fr_unc_off fr_unc_sig N_crlb_max N_crlb_min N_tol_max N_tol_min off_crlb_max off_crlb_min off_max off_min sig_crlb_max sig_crlb_min index fr_N_t fr_off_t fr_sigx_t fr_sigy_t N_temp N_crlb_temp xf_all_temp xf_crlb_temp yf_all_temp yf_crlb_temp off_all_temp off_crlb_temp llv_temp lp_temp framenum_all_temp sigx_temp sigy_temp sigy_crlb_temp sigx_crlb_temp

%% save workspace and append '_tol' to file name
save([pwd,'\Tol\',bname(1:end-4), '_tol']);