% Plot yer Regions
%
% This script will plot the regions of a given index from a selected
% toleranced file from the image_analysis_spt function
%  This will also create a convex hull to show you the "shape" of the
%  region identified by the points
%

clearvars
close all;
clc;

index = 10;


[fname, fpath] = uigetfile('*tol.mat');

load([fpath,fname]);
figure
for i = 1:max(number)
    ind = find(number == i);
    
    try
    k = convhull(xf_all(ind)*q,yf_all(ind)*q);
    xcm = sum(xf_all(ind)*q)/numel(ind);
    ycm = sum(yf_all(ind)*q)/numel(ind);
    
    
    subplot(2,5,i);
    hold on
    plot(xf_all(ind(k))*q,yf_all(ind(k))*q);
    plot(xf_all(ind)*q,yf_all(ind)*q,'.b');
    plot(xcm,ycm,'.r')
    title(['Number ', num2str(i)]);
    hold off
    catch lsterr
    end
end