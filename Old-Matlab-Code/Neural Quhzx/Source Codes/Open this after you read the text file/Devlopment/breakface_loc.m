%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Breakface Localization
%
% A new localization approach based around brute force localization
%
% AJN 7/17/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clearvars; close all;

[fname, fpath] = uigetfile('*.tif');

% Load Image
i1 = readtiff([fpath,fname]);

% Background Subtract
iprod = single(rollingball(double(i1)));
[m,n,o] = size(iprod);
xf_all = [];
yf_all = [];
N = [];
sigx_all = [];
sigy_all = [];
off_all = [];
xf_crlb = [];
yf_crlb = [];
N_crlb = [];
sigx_crlb = [];
sigy_crlb = [];
off_crlb = [];
llv = [];
framenum_all = [];

% Localize all pixels
for i = 1:o
    i/o
    [xf,xc, yf,yc, Np,  Nc, sigx, sigxc, sigy, sigyc,off, offc, lv] = break_loc(iprod(:,:,i));
       xf_all = [   xf_all; xf];
       yf_all = [   yf_all; yf];
            N = [        N; Np];
     sigx_all = [ sigx_all; sigx];
     sigy_all = [ sigy_all; sigy];
      off_all = [  off_all; off];
      xf_crlb = [  xf_crlb; xc];
      yf_crlb = [  yf_crlb; yc];
       N_crlb = [   N_crlb; Nc];
    sigx_crlb = [sigx_crlb; sigxc];
    sigy_crlb = [sigy_crlb; sigyc];
     off_crlb = [ off_crlb; offc];
          llv = [      llv; lv];
 framenum_all = [framenum_all ; i*ones(numel(xf),1)];
end

save([fname(1:end-4),'_broke.mat'],'xf_all', 'yf_all', 'xf_crlb', 'yf_crlb', 'N', 'N_crlb', 'off_all', 'off_crlb', 'sigx_all', 'framenum_all', 'sigx_crlb', 'sigy_all', 'sigy_crlb', 'llv');