function [fits] = wave_fit(i1)
% This is a script to find the 'Fitting Parameters' of wavelet filtered
% images, initially this calculates CoM for all images, but will eventually
% include the 3D correlation metric.
% for an image of mxmxo the output will be a ox3 array [xcm, ycm, zcr]
% AJN 9/20/18

[m,n,o] = size(i1);
wind = (1:m) - (m+1)/2;
[X,Y] = meshgrid(wind);
fits = zeros(o,4);
for i = 1:o % perform 'fit' over each image
    si1 = sum(sum(i1(:,:,i)));
    fits(i,1) = sum(sum(X.*i1(:,:,i)/si1));
    fits(i,2) = sum(sum(Y.*i1(:,:,i)/si1));
    fits(i,3) = si1;
end