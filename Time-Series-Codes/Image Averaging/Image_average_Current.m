% Image Averager
%
% A simple script to take in several time-averaged trials and create 1
% final averaged copy
%
% Ajn 3/7/17
%
% Version 2 will automatically combine all files in a folder and
% additionally output a single image of the peak frame - the average frames
clear all; close all; clc;

baseline = 100;
sp = 110;
psize = 3;

[fname, fpath] = uigetfile('*.tif','Choose a baseline tiff'); % Allow user to select a director with images in it

cd(fpath); % change working directory to chosen directory
files = dir('*.tif'); % get a list of all images in directory
indstr = 'phluorin';
for i = 1:numel(files)
    lstri = files(i).name;
    if contains(lstri,indstr)
        bigind = i;
        break
    end
end
flag = 0;
% clear bigind

if exist('bigind')
    [peakresponse] = getpeakresponse(files(bigind).name, psize);
    flag = 1;
end
clear i1
x = [];
if exist(['measured_results0.mat'])
    sn = input('Use previously measured positions?','s');
else
    sn = 'n';
end
for j = 1:numel(files) % loop over every image found
    finfo = dir([files(j).name(1:end-4),'*']); % find images with the same base file name, as andor labels _1, _2, ...
    % for additional items with
    % the same name this will
    % behave conveniently
    if numel(finfo) > 1  % only perform averaging if there are more than 1 file with the base name
        for i = 1:numel(finfo)
            %% Image Preparation
            im1 = readtiff(finfo(i).name);
%             iprod = rollingball(im1,4,4);
            
%             clear im1
%             base = mean(iprod(:,:,1:baseline-1),3);
            [iprod] = fsub(im1,sp);
            if flag == 0
                peakresponse = 1;
            end
            base = mean(iprod(:,:,1:baseline-1),3);
            % iprod - base gives deltaf
            func_time_series_2(iprod-base, peakresponse, psize, i, sn);
       end
    end
end

clear fluor
avef = disp_data(1.2);