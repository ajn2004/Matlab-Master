% Pix 2 photon measurement
%
% This will take a tiff stack and return the graph of mean / variance of
% every pixel, then fit a line to the data, the slope of the line being the
% pixel to photon ratio
%
% AJN 3/15/17
clear all; close all; clc;
% warning('off');
% File selection and loading
[fnm, fpath] = uigetfile('*tif');
cd(fpath);
files = dir('*.tif');
% imagfo = imfinfo([fpath, fname]);
vars = [];
aves = [];
for i = 1:numel(files)
i1 = readtiff(files(i).name);
%  imag = fitsinfo(fnm); % get image about file
%     i1 = fitsread(fnm,'Info', imag);

vars = cat(3,vars,var(i1,0,3));
aves = cat(3,aves,mean(i1,3));
end
vars = mean(vars,3);
aves = mean(aves,3);
lim = 12000;
fits = polyfit(aves(aves<lim),vars(aves<lim),1);

disp(['The pixel to photon ratio is ', num2str(fits(1))]);

plot(aves(aves<lim),vars(aves<lim),'.b')
hold on
plot([min(aves(aves<lim)), max(aves(aves<lim))],[fits(1)*min(aves(aves<lim))+fits(2), fits(1)*max(aves(aves<lim))+fits(2)],'r','LineWidth', 2);
hold off
xlabel('Mean value')
ylabel('Variance');