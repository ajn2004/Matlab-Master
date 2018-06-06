%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation to test pix2photon_conversion script
% AJN 7/4/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace
clear all
close all
clc


%% Parameters of map of interest
size = 50;  %size in square pixels
depth = 5000;
im_num = 10; %number of files to be created

%% Pixel map of pix to photo ratios
% pi2ph = 35+rand(size,size)*3;           %builds pixel map
pi2ph = 35+20*rand(size,size);

%% Builds file stacks
for k = 1:im_num
    clear A
    A = poissrnd(1+k*20,size,size,depth);   %creates array of poisson distributed values
    for l = 1:depth
        Apix(:,:,l)= A(:,:,l).*pi2ph;      %converts each frame from photons to pixels
    end
    
%% Writes a tiff stack

outputFileName = strcat('tiff_sim_',num2str(k),'.tif');
for K=1:numel(Apix(1, 1, :))
    
   imwrite(uint16(Apix(:,:,K)), outputFileName, 'tif','WriteMode', 'append');
end
end
