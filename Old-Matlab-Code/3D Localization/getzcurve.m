%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Z-Curve
%
% This program will extract sigx / sigy information and attempt to
% ascertain a "height" from a naming convetion. This program is designed to
% run on a 3D calibration set where each file is a series of images at a
% specific position
% AJN 1-9-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
clc;

files = dir('*.mat');


sigxa = [];
sigxc = [];
sigyc = [];
sigya = [];
height = [];
iloca =[];
sigx = [];
sigy = [];
ssigx =[];
ssigy = [];
heights = [];
xf = [];
yf = [];
zf = [];
for i = 1:numel(files)
    load(files(i).name);
    xf = [xf;xf_all];
    yf = [yf;yf_all];
%     [xf_all,xf_crlb, yf_all,yf_crlb, N, N_crlb, sigx_all, sigx_crlb,sigy_all, sigy_crlb,off_all, off_crlb, framenum_all, llv, y] = da_locs(sum(iloc,2), 1, [0, 0], 0);
    sigxa = [sigxa;sigx_all];
    sigx = [sigx;mean(sigx_all)];
    ssigx = [ssigx;std(sigx_all)/numel(sigx_all)^0.5];
%     sigxc = [sigxc;rx0_crlb];
    sigya = [sigya;sigy_all];
    sigy = [sigy;mean(sigy_all)];
    ssigy = [ssigy;std(sigy_all)/numel(sigx_all)^0.5];
%     sigyc = [sigyc;sigy_crlb];
%     iloca = [iloca,iloc(:,ind)];
    inds = strfind(files(i).name,'_');
    heights  =[heights; str2num(files(i).name(inds(1)+1:inds(2)-1))];
    height = [height; ones(numel(xf_all),1)*str2num(files(i).name(inds(1)+1:inds(2)-1))];
    clearvars -except files sigxa ssigx ssigy sigya iloca fnuma i sigxc sigyc height heights sigx sigy xf yf
end
z = 50*(heights-19);
errorbar(heights*50,sigx,ssigx,'.');
hold on
errorbar(heights*50,sigy,ssigy,'.');
hold off
legend('Sigx','Sigy');
