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

baseline = 50;

[fname, fpath] = uigetfile('*.tif'); % Allow user to select a director with images in it

cd(fpath); % change working directory to chosen directory
files = dir('*.tif'); % get a list of all images in directory
for j = 1:numel(files) % loop over every image found
    finfo = dir([files(j).name(1:end-4),'*']); % find images with the same base file name, as andor labels _1, _2, ... 
                                               % for additional items with
                                               % the same name this will
                                               % behave conveniently
    if numel(finfo) > 1  % only perform averaging if there are more than 1 file with the base name
        %% Image Averaging
        imagf1 = imfinfo(files(j).name); % get image information of the base file, which should be identical over all subsequent images in the set
        i2 = zeros(imagf1(1).Height,imagf1(1).Width,numel(imagf1)); % preallocate space
        for i = 1:numel(finfo)  % Loop over all images found
            i
            imagfo = imfinfo(finfo(i).name); % get image information of current image
            for k = 1:numel(imagfo) % Loop over all frames in image
                i1(:,:,k) = double(imread(finfo(i).name,k)); % load image into memory
            end
            i2 = i2 +i1; % add image to preallocated space
        end
        i2 = (i2./numel(finfo)); % divide final total image by number of images used (i.e. average)
        ave2 = mean(i2(:,:,1:baseline-1),3);
        ifin = round(i2);
        imwrite(uint16(ifin(:,:,1)),['Averaged_',files(j).name]); % start writing averaged tiff
        count = 2; % counting variable to get around tiff writing bug in matlab
        
        %% Tiff Writing Areas
        while true % Continuously repeat until broken
            try % Attempt to write image, if unsuccesful this will go to the catch and the loop will repeat until frame is written
                imwrite(uint16(ifin(:,:,count)),['Averaged_',files(j).name],'writemode','append');
                count = count +1; %counting variable updated after successful frame writing
                if count > numel(imagf1) % if counting variable gets to frames + 1 break the while statement
                    break
                end
            catch lsterr
            end
        end
        
        %% Peak Image Section
        % The logic of this section is that the brightest frame should have
        % the brightest pixels over the whole frame, so by summing all
        % pixels on each frame independtly we should be able to find the
        % peak response frame by looking for the largest value
        
        sums = sum(sum(i2,1),2); % sum all frames and store the final result
        peak = find(sums == max(sums)); % find the largest value index which will result in the largest frame
        
        avefin = mean(i2,3);  % take the average over all frames
        peakfin = i2(:,:,peak) - avefin; % Create the image of differences
        imwrite(uint16(peakfin), ['Sub_',files(j).name]); % Write that to the folder
        
    end
end