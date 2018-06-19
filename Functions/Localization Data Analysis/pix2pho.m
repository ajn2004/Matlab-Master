% Pix 2 photon measurement
%
% This will take a tiff stack and return the graph of mean / variance of
% every pixel, then fit a line to the data, the slope of the line being the
% pixel to photon ratio
%
% AJN 3/15/17
clear all; close all; clc;
warning('off');
% File selection and loading
[fnm, fpath] = uigetfile('*tif');
cd(fpath);
% imagfo = imfinfo([fpath, fname]);
% i1 = readtiff([fpath,fname]);
 imag = fitsinfo(fnm); % get image about file
    i1 = fitsread(fnm,'Info', imag);
lim = 12000;
vars = var(i1,0,3);
aves = mean(i1,3);

fits = polyfit(aves(aves<lim),vars(aves<lim),1);

disp(['The pixel to photon ratio is ', num2str(fits(1))]);

plot(aves(aves<lim),vars(aves<lim),'.b')
hold on
plot([min(aves(aves<lim)), max(aves(aves<lim))],[fits(1)*min(aves(aves<lim))+fits(2), fits(1)*max(aves(aves<lim))+fits(2)],'r','LineWidth', 2);
hold off
xlabel('Mean value')
ylabel('Variance');