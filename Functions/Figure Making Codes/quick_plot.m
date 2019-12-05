clearvars; clc; close all;
% [fname, fpath] = uigetfile('*dast*');
files = dir('*dast*');
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\Channel_net');
for i = 1:numel(files)
    fname = files(i).name;
load(fname);
try
zf_all = func_shift_correct(ncoords(:,3),framenumber,1);
xf_all = ncoords(:,1);
yf_all = ncoords(:,2);
ID = xf_all < 175;
ind = logical(1-ID);

xf_blue = xf_all(ID);
yf_blue = yf_all(ID);
xf_red = xf_all(ind);
yf_red = yf_all(ind);
zf_blue = zf_all(ID);
zf_red = zf_all(ind);
p_coords = [xf_blue,yf_blue];
xf_new = xnet(p_coords.');
yf_new = ynet(p_coords.');
xf_fin = xf_new(:);
yf_fin = yf_new(:);

plot(xf_fin*q,q*yf_fin,'.b');
hold on
plot(xf_red*q,q*yf_red,'.r')
legend('Cono Halo','Munc13-Dendra2');
xlabel('Microns')
ylabel('Microns')
title('Cono-Halo and Munc13-Dendra2')
hold off
waitforbuttonpress
catch lsterr
end
end