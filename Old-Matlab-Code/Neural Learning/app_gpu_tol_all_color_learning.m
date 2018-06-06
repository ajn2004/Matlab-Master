function [x, y] = app_gpu_tol_all_color_learning(it, dir_name, base_name, sig_min, sig_max, sig_crlb_min, sig_crlb_max, fr_unc_sig,  N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off)
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
y = [];
x = [];
try
    N_temp = N;
    N_crlb_temp = N_crlb;
    xf_all_temp = xf_all;
    xf_crlb_temp = xf_crlb;
    yf_all_temp = yf_all;
    yf_crlb_temp = yf_crlb;
    off_all_temp = off_all;
    off_crlb_temp = off_crlb;
    sigx_all_temp = sigx_all;
    sigy_all_temp = sigy_all;
    sigx_crlb_temp = sigx_crlb;
    sigy_crlb_temp = sigy_crlb;
    llv_temp = llv;
    framenum_all_temp = framenum_all;
    fr_N = N_crlb.^0.5./N;
    fr_off = off_crlb.^0.5./off_all;
    fr_sigx = sigx_crlb.^0.5./sigx_all;
    fr_sigy = sigy_crlb.^0.5./sigy_all;
    
    %% load multicolor info if it exists
    if exist('nrat');
        nrat_temp = nrat;
        ratio_temp = ratio;
        clear nrat ratio
    end
    %% clear variables
    clear N N_crlb xf_all xf_crlb yf_all yf_crlb off_all off_crlb llv lp framenum_all sigx_all sigy_all sigx_crlb sigy_crlb
    %% determine index of acceptable molecules
%     if it > 1
    index = find(sigx_crlb_temp > sig_crlb_min & sigy_crlb_temp > sig_crlb_min & sigx_crlb_temp < sig_crlb_max & sigy_crlb_temp < sig_crlb_max & sigx_all_temp < sig_max & sigy_all_temp < sig_max & sigx_all_temp > sig_min & sigy_all_temp > sig_min & fr_sigx <fr_unc_sig & fr_sigy <fr_unc_sig & round(xf_all_temp) <= xc + 1 & round(xf_all_temp) >= xc - 1 & round(yf_all_temp) <= yc + 1 & round(yf_all_temp) >= yc -1 & fr_N < fr_unc_N & fr_off < fr_unc_off & N_temp < N_tol_max & N_temp > N_tol_min & N_crlb_temp > N_crlb_min & N_crlb_temp < N_crlb_max & xf_crlb_temp < xf_crlb_max & xf_crlb_temp > xf_crlb_min & yf_crlb_temp < yf_crlb_max & yf_crlb_temp > yf_crlb_min & off_crlb_temp > off_crlb_min & off_crlb_temp < off_crlb_max & off_all_temp >off_min & off_all_temp <off_max);
    index0 =find(sigx_crlb_temp < sig_crlb_min | sigy_crlb_temp < sig_crlb_min | sigx_crlb_temp > sig_crlb_max | sigy_crlb_temp > sig_crlb_max | sigx_all_temp > sig_max | sigy_all_temp > sig_max | sigx_all_temp < sig_min | sigy_all_temp < sig_min | fr_sigx >fr_unc_sig | fr_sigy >fr_unc_sig | round(xf_all_temp) >  xc + 1 | round(xf_all_temp) <  xc - 1 | round(yf_all_temp)  > yc + 1 | round(yf_all_temp) <  yc - 1 | fr_N > fr_unc_N | fr_off > fr_unc_off | N_temp > N_tol_max | N_temp < N_tol_min | N_crlb_temp < N_crlb_min | N_crlb_temp > N_crlb_max | xf_crlb_temp > xf_crlb_max | xf_crlb_temp < xf_crlb_min | yf_crlb_temp > yf_crlb_max | yf_crlb_temp < yf_crlb_min | off_crlb_temp < off_crlb_min | off_crlb_temp > off_crlb_max | off_all_temp <off_min | off_all_temp >off_max);
% %     else
%     index = find(sigx_crlb_temp > sig_crlb_min & sigy_crlb_temp > sig_crlb_min & sigx_crlb_temp < sig_crlb_max & sigy_crlb_temp < sig_crlb_max & sigx_all_temp < sig_max & sigy_all_temp < sig_max & sigx_all_temp > sig_min & sigy_all_temp > sig_min & fr_sigx <fr_unc_sig & fr_sigy <fr_unc_sig & fr_N < fr_unc_N & fr_off < fr_unc_off & N_temp < N_tol_max & N_temp > N_tol_min & N_crlb_temp > N_crlb_min & N_crlb_temp < N_crlb_max & xf_crlb_temp < xf_crlb_max & xf_crlb_temp > xf_crlb_min & yf_crlb_temp < yf_crlb_max & yf_crlb_temp > yf_crlb_min & off_crlb_temp > off_crlb_min & off_crlb_temp < off_crlb_max & off_all_temp >off_min & off_all_temp <off_max);
%     index0 =find(sigx_crlb_temp < sig_crlb_min | sigy_crlb_temp < sig_crlb_min | sigx_crlb_temp > sig_crlb_max | sigy_crlb_temp > sig_crlb_max | sigx_all_temp > sig_max | sigy_all_temp > sig_max | sigx_all_temp < sig_min | sigy_all_temp < sig_min | fr_sigx >fr_unc_sig | fr_sigy >fr_unc_sig | fr_N > fr_unc_N | fr_off > fr_unc_off | N_temp > N_tol_max | N_temp < N_tol_min | N_crlb_temp < N_crlb_min | N_crlb_temp > N_crlb_max | xf_crlb_temp > xf_crlb_max | xf_crlb_temp < xf_crlb_min | yf_crlb_temp > yf_crlb_max | yf_crlb_temp < yf_crlb_min | off_crlb_temp < off_crlb_min | off_crlb_temp > off_crlb_max | off_all_temp <off_min | off_all_temp >off_max);
% %     end
    %% redefine variables to be only of acceptable molecules
    y = ones(numel(index),1);
    y = [y;zeros(numel(index0),1)];
catch lsterr
end




%% build x image matrix

try
    % positive images
    disp('Building Positive Images');
    for i = 1:numel(index)
        i1 = i5(:,:,index(i));
        x = [x;i1(:).'];
    end
    
    % negative passed loc images
    disp('Building Negative Images');
    for i = 1:numel(index0)
        i1 = i5(:,:,index0(i));
        x = [x;i1(:).'];
    end
    
    % add in i6;
catch lsterr
end
disp('Building Final Images');
if numel(i6(1,1,:)) > 5000
    y = [y;zeros(5000,1)];
    indy = randperm(numel(i6(1,1,:)));
    for i = 1:5000
        i1 = i6(:,:,indy(i));
        x = [x; i1(:).'];
    end
else
    y = [y;zeros(numel(i6(1,1,:)),1)];
    for i = 1:numel(i6(1,1,:))
        i1 = i6(:,:,indy(i));
        x = [x; i1(:).'];
    end
end