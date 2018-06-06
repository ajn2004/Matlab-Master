%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Loraxian Viewer
%
% This is a script for viewing the hidden layers of the loraxian transform
% neural network
%
% Andrew Nelson
%
%12/22/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

[fname, fpath] = uigetfile('*.mat');
load([fpath,fname]);

%initialize final image space as 2 5x5s
i1 = zeros(63)-1;
i2 = i1;

for i = 1:numel(theta1(:,1))/2
    im1 = reshape(theta1(i,2:end),9,9);
    im1s = im1./max(im1);
    im1f(:,:,i) = im1s - min(im1s);
end

for i = numel(theta1(:,1))/2+1:numel(theta1(:,1))
    im1 = reshape(theta1(i,2:end),9,9);
    im1s = im1./max(im1);
    im2f(:,:,i-numel(theta1(:,1))/2) = im1s - min(im1s);
end

for i =1:5
    for j = 1:5
        i1(4+(i-1)*12:i*12,4+(j-1)*12:j*12) = im1f(:,:,(i-1)*5 + j);
        i2(4+(i-1)*12:i*12,4+(j-1)*12:j*12) = im2f(:,:,(i-1)*5 + j);
    end
end

subplot(1,2,1);imagesc(i1);colormap('gray')

subplot(1,2,2);imagesc(i2);colormap('gray');