%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Da Storm
%
% The newest localization algorithm based off quhzx. Instead of neural
% network identification we will be using image PCA to identify the
% important areas
%
% AJN 7/24/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clearvars
clc
p = mfilename('fullpath');
[fpath, fname, fext] = fileparts([p, '.m']);
addpath(fpath);
addpath([fpath,'\da_c']);
addpath([fpath,'\da_functions']);
pixw = 3;
pix2pho = 33.4;
q = 0.127;
% Grab File
[fname, data_d] = uigetfile('*.tif');
cd(data_d)

% Convert Variabls
pix2pho = single(pix2pho);
q = single(q);
% Load file
i1 = readtiff([data_d, fname])/pix2pho;
[m,n,o] = size(i1);
% Rolling Ball Background Subtract
iprod = rollingball(i1);

% Peak Detection
ip = ImsPCA(iprod,3);
[thresh] = find_thresh(ip);
% varys = {'xf_a', 'xf_c','yf_a', 'yf_c', 'N_a','N_c','off_a','off_c','sigx_a','sigx_c','sigy_a','sigy_c','llv','fnum_a'};
% varys2 = {'xf_o', 'xf_co','yf_o', 'yf_co', 'N_o','N_co','off_o','off_co','sigx_o','sigx_co','sigy_o','sigy_co','llv_o'};
% for i = 1:numel(varys)
%     eval([varys{i},'=[];']);
% end

dps = da_peaks(ip,thresh);
clear ip


% divide up the data
[iloc, fnum, cents] = divide_up(iprod, pixw, dps);

[xf_o,xf_co, yf_o,yf_co, N_o, N_co, sigx_o, sigx_co,sigy_o, sigy_co,off_o, off_co, fnout, llv_o] = da_locs(iloc, fnum, cents);
save([fname(1:end-4),'_dast.mat'], 'xf_o' , 'xf_co' , 'yf_o' , 'yf_co' , 'N_o' , 'N_co' , 'sigx_o' , 'sigx_co' , 'sigy_o' , 'sigy_co' ,'off_o' , 'off_co', 'fnout', 'llv_o','thresh','pixw','q','pix2pho');
