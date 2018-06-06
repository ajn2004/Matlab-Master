%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show Tol for GPU
%
% Shows the result of applying specific tolerances for a GPU analyzed file
% This code has gusto!
% AJN 2/18/16
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% clear all
% clc


if ~ exist('xf_all')

    [fname, fpath] = uigetfile('*.mat');
    load([fpath,fname]);
    cd(fpath)
end
% Tolerance Variables
     iln_min = -1;
   N_tol_min = 250;        % minimum number of photons
   N_tol_max = 700000;       % maximum number of photons
     off_min = 0;          % minimum number of offset photons
     off_max = 100;         % maximum number of offset photons
  N_crlb_min = 0;          % minimum variance in number of photons
  N_crlb_max = 100000;      % maximum variance in number of photons
     max_unc = 30;       % maximum error in number of x-position
%     max_zunc = 110;      % max reported uncertainty in z
     sig_min = 100;     % sigma minimum in nm
     sig_max = 500;     % sigma maximum in nm
     sig_crlb_min = 0;
     sig_crlb_max = 10;
off_crlb_min = 0.0;          % minimum error in number of offset photons
off_crlb_max = 5;          % maximum error in number of offset photons
      zf_max = 600;
      zf_min = -500;
    fr_unc_N = 0.125;        % fractional uncertainty in N
  fr_unc_off = 0.4;        % Fractional uncertainty in offset
      fr_sig = 0.4;         % fractional uncertainty for sigmas

zf_all = getdz(sigx_all,sigy_all)/(q);

% inloc = iloc(:,logical(y));
%% Tolerancing


xf_crlb_max = (max_unc/(q*1000))^2;
yf_crlb_max = (max_unc/(q*1000))^2;
% zf_crlb_max = (max_zunc/(q*1000))^2;
iln = llv./N;
fr_N = N_crlb.^0.5./N;
fr_off = off_crlb.^0.5./off_all;
fr_sx = sigx_crlb.^0.5./sigx_all;
fr_sy = sigy_crlb.^0.5./sigy_all;
N_temp = N;
N_crlb_temp = N_crlb;
xf_all_temp = xf_all;
xf_crlb_temp = xf_crlb;
yf_all_temp = yf_all;
yf_crlb_temp = yf_crlb;
zf_all_temp = zf_all;
% zf_crlb_temp = zf_crlb;
sigx_all_temp = sigx_all;
sigx_crlb_temp = sigx_crlb;
sigy_all_temp = sigy_all;
sigy_crlb_temp = sigy_crlb;
off_all_temp = off_all;
off_crlb_temp = off_crlb;
llv_temp = llv;
framenum_all_temp = framenum_all;

zf_max = zf_max/(q*1000);
zf_min = zf_min/(q*1000);
sig_min = sig_min/(q*1000);
sig_max = sig_max/(q*1000);

% index = find(xf_all >= 5 & yf_all >= 5 & xf_all <= max(xf_all) - 5 &...
%     yf_all <= max(yf_all) - 5 & iln >= iln_min & (zf_all) <= zf_max & zf_all >= zf_min &...
%     fr_N <= fr_unc_N & fr_off <= fr_unc_off & N_temp <= N_tol_max &...
%     N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max &...
%     xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max &...
%     off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max &...
%     off_all_temp >= off_min & off_all_temp <= off_max & sigx_all_temp <= sig_max...
%     & sigx_all_temp >= sig_min & sigy_all_temp <= sig_max & sigy_all_temp >= sig_min...
%     & sigx_crlb_temp <= sig_crlb_max & sigx_crlb_temp >= sig_crlb_min...
%     & sigy_crlb_temp <= sig_crlb_max & sigy_crlb_temp >= sig_crlb_min...
%     & fr_sx <= fr_sig & fr_sy <=fr_sig);
index = find(xf_all >= 5 & yf_all >= 5 &...
    iln >= iln_min & (zf_all) <= zf_max & zf_all >= zf_min &...
    fr_N <= fr_unc_N & fr_off <= fr_unc_off & N_temp <= N_tol_max &...
    N_temp >= N_tol_min & N_crlb_temp >= N_crlb_min & N_crlb_temp <= N_crlb_max &...
    xf_crlb_temp <= xf_crlb_max & yf_crlb_temp <= yf_crlb_max &...
    off_crlb_temp >= off_crlb_min & off_crlb_temp <= off_crlb_max &...
    off_all_temp >= off_min & off_all_temp <= off_max & sigx_all_temp <= sig_max...
    & sigx_all_temp >= sig_min & sigy_all_temp <= sig_max & sigy_all_temp >= sig_min...
    & sigx_crlb_temp <= sig_crlb_max & sigx_crlb_temp >= sig_crlb_min...
    & sigy_crlb_temp <= sig_crlb_max & sigy_crlb_temp >= sig_crlb_min...
    & fr_sx <= fr_sig & fr_sy <=fr_sig);
index0 =find( xf_all < 2 | yf_all < 2 | xf_all > max(xf_all) - 2 |...
    yf_all > max(yf_all) - 2 | iln < -1 | abs(zf_all) >  zf_max |  fr_N > fr_unc_N |...
    fr_off > fr_unc_off | N_temp > N_tol_max | N_temp < N_tol_min | N_crlb_temp < N_crlb_min |...
    N_crlb_temp > N_crlb_max | xf_crlb_temp > xf_crlb_max |...
    yf_crlb_temp > yf_crlb_max | off_crlb_temp < off_crlb_min |...
    off_crlb_temp > off_crlb_max | off_all_temp <off_min | off_all_temp >off_max);
% index = find(llv./N >= -1);
% index0 = find(llv./N <-1);


figure('units','normalized','outerposition',[0.5 0 0.5 1])
subplot(3,2,1); hist(N(index), 2*numel(N(index))^(1/3)); title('Histogram of N');
subplot(3,2,3); hist(fr_N(index), 2*numel(N(index))^(1/3)); title('Histogram of fr N');
subplot(3,2,2); hist(off_all(index), 2*numel(N(index))^(1/3)); title('Histogram of Offset');
subplot(3,2,4); hist(fr_off(index), 2*numel(N(index))^(1/3)); title('Histogram of fr off');
subplot(3,2,5); hist(N_crlb(index), 2*numel(N(index))^(1/3)); title('Histogram of N crlb');
subplot(3,2,6); hist(off_crlb(index), 2*numel(N(index))^(1/3)); title('Histogram of off crlb');
subplot(4,2,8); hist(sigx_all(index)*q*1000*2, 2*numel(N(index))^(1/3)); title('Histogram of sigma x ');xlabel('Radius in nm'); ylabel('Frequency');
subplot(4,2,7); hist((sigy_all(index).*sigx_all(index)).^0.5*2*q*1000, 2*numel(N(index))^(1/3)); title('Histogram of r0 all*2');xlabel('Radius in nm'); ylabel('Frequency');


figure('units','normalized','outerposition',[0 0.5 0.5 0.5])
subplot(1,2,1);
hold on
plot(xf_all(index0),yf_all(index0),'.r');
plot(xf_all(index),yf_all(index),'.b');
hold off
legend('Exclude','Include')
title(['Number of molecules in tolerance =', num2str(numel(xf_all(index)))])
xlabel([num2str(numel(xf_all(index))/numel(xf_all)),'% of total localizations']);
axis image
subplot(1,2,2);
% imagesc(reshape(mean(inloc(:,index),2),11,11));
histogram(zf_all(index)*q*1000,'Normalization','Probability');
hold on
histogram(zf_all(index0)*q*1000,'Normalization','Probability');
hold off
legend('inside','outside');
title('Histogram of included axial positions');
figure('units','normalized','outerposition',[0 0 0.5 0.5])
subplot(2,2,1); hist(xf_crlb(index).^0.5*q*1000, 2*numel(N(index))^(1/3)); title('Uncertainty in X');xlabel('Uncertainty in nm'); ylabel('Frequency');
% subplot(2,2,3); hist(zf_crlb(index).^0.5*q*1000, 2*numel(N(index))^(1/3)); title('Uncertainty in Z');xlabel('Uncertainty in nm'); ylabel('Frequency');
subplot(2,2,2); hist(yf_crlb(index).^0.5*q*1000, 2*numel(N(index))^(1/3)); title('Uncertainty in Y');xlabel('Uncertainty in nm'); ylabel('Frequency');
subplot(2,2,3); histogram(q*1000*sigx_all(index)); hold on; histogram(q*1000*sigy_all(index)); title('Sigma Values'); legend('Sig x', 'Sig y'); xlabel('Width in nm');
subplot(2,2,4); histogram(sigx_crlb(index).^0.5./sigx_all(index)); hold on; histogram(sigy_crlb(index).^0.5./sigy_all(index)); title('Sigma frac unc'); legend('Sig x', 'Sig y'); xlabel('Width in nm');
whitebg([1 1 1])

save('C:\Users\AJN Lab\Desktop\Code Development\Matlab Testing Folder\Neural Quhzx\GPU Tolerance\these_tol.mat', 'zf_min','N_tol_min','N_tol_max','off_min','off_max','N_crlb_min','N_crlb_max','sig_crlb_min','sig_max','sig_min','fr_sig','sig_crlb_max','xf_crlb_max','yf_crlb_max','off_crlb_min','off_crlb_max','iln_min','fr_unc_N','fr_unc_off','zf_max');