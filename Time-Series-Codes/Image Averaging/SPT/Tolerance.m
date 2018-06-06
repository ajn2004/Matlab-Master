%% User Specified Tolerances
% choose the min/max values based on the histograms of these variables
% i.e. N, off_all, xf_crlb, yf_crlb, N_crlb, off_crlb, lp, and llv
% Tolerance Variables
N_tol_min = 50;         % minimum number of photons
N_tol_max = 50000;       % maximum number of photons
off_min = 0;            % minimum number of offset photons
off_max = 100;           % maximum number of offset photons
N_crlb_min = 0;         % minimum variance in number of photons
N_crlb_max = 20000;       % maximum variance in number of photons
sigma2_min = 150;        % minimum sigma value
sigma2_max = 750;         % maximum sigma value
xf_crlb_min = 0;        % minimum error in number of x-position
xf_crlb_max = 0.4;      % maximum error in number of x-position
yf_crlb_min = 0;        % minimum error in number of y-position
yf_crlb_max = 0.4;      % maximum error in number of y-position
off_crlb_min = 0;     % minimum error in number of offset photons
off_crlb_max = 50;     % maximum error in number of offset photons
sig_crlb_min = 0;
sig_crlb_max = 2;
llv_max = -1;
llv_min = -100;
fr_unc_N = 0.8;           % fractional uncertainty in N
fr_unc_off = 1;       % Fractional uncertainty in offset
fr_unc_sig = 0.8;       % Fractional uncertainty in width
[fname, fpath] = uigetfile('*results.mat');
q = 0.156;
cd(fpath);

finfo = dir('*results.mat');

for i = 1:numel(finfo)
%     load(finfo(i).name);
    func_tol_spt(finfo(i).name,q, llv_max, llv_min, sigma2_min, sigma2_max, sig_crlb_min, sig_crlb_max, fr_unc_sig, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off);
end
