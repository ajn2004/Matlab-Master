%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative Selector
%
% The neural net used in image_neural_2.cu is effective at picking out
% molecules, but with that it picks up a lot of noise. This program allows
% a user to select regions of an image to be used in a negative example
% training set. 
%
% AJN 12/17/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

numex = 10;

%% Basic startup and file selection
[fname, fpath] = uigetfile('*.tif', 'Select an image to choose from');
addpath = pwd;
boxw = 3;
cd(fpath);
numfd = 0;
imagfo = imfinfo([fpath,fname]);
xf = [];
yf = [];
frames = [];
% continue to look over frames until number of examples is reached
while numfd < numex
    frame = randi(numel(imagfo));
    i1 = imread([fpath,fname], frame, 'Info',imagfo);
    imagesc(i1, [min(i1(:)), max(i1(:))]);
    axis equal
    colormap('gray');
    title('Select a pixel that does not have a molecule in it, press enter when done');
    [x, y] = ginput(numex - numfd);
    for i = 1:numel(x);
        numfd = numfd + 1;
        xind = ceil(x(i))-3:ceil(x(i))+3;
        yind = ceil(y(i))-3:ceil(y(i))+3;
        i2(:,:,numfd) = i1(yind,xind);
    end
    frames = [frames; frame*ones(numel(x),1)];
    x = ceil(x);
    y = ceil(y);
end

%% at this point we have a series of regions
    
        
            
        