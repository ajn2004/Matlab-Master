function [pixmap] = pix2photo_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to measure the pixel to photon ratio of the sCMOS camera
% 7/7/13 AJN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE TO USER: This program assumes that you have taken several images of
% a homegenous light source at different intensities and have labeled them
% sequentially such as example_1.tif. 
% This program will automatically load all files in sequence and analyze them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the workspace
clear all
close all
clc

%% Declare some variables
file_num = input('How many files will be used to analyze the pixel to photon map? ');  % number of files in your series
% frame = 5000;    % number of frames in each picture

%% Load primary file and creat series setup
[fname fpath] = uigetfile('*.tif', ' Choose the file that starts the series');

sname = fname(1:numel(fname)-5);        %creates a name for all series

frame = numel(imfinfo(strcat(fpath,fname)));
%% Creation of mean pixel map and average pixel map
for k = 0:file_num-1
    
    %% Load all files in series
    for i = 1:frame
        A(:,:,k*frame+i) = imread(strcat(fpath,sname,num2str(k+1),'.tif'),i);
    end
    %% Creation of average and variance of each pixel
    Aave(:,:,k+1) = mean(double(A(:,:,k*frame+1:(k+1)*frame)),3);       % average
    Avar(:,:,k+1) = var(double(A(:,:,k*frame+1:(k+1)*frame)),0,3);      % variance

end


%% Calculation of pixel to photon ratio
for m = 1:length(A(:,1,1))
    for n = 1:length(A(1,:,1))
        clear p
        p = polyfit(Aave(m,n,:),Avar(m,n,:),1); % uses polyfit to measure slope of mean vs. variance
        pixmap(m,n) = p(1);                     % generates a pixel map of slopes (pix2pho ratio)
    end
end

save(strcat(fpath,'pix2photon_map'), pixmap)

%% For simulation purposes
% mean(mean(abs(pixmap-pi2ph)./pi2ph))