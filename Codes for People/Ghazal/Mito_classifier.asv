clearvars; close all; clc; % clean up

bkgn_p = 20;  % percentage to use in background calculation

[mname, mpath] = uigetfile('*mtRFP*');  % Selecting file in desired folder
cd(mpath);  % change working directory to selected folder
clear mname mpath % clear file selection variables
files = dir('*mtRFP*'); % get file list of all mitochondria images

totrat = []; % preallocate total ratio variable
for l = 1:numel(files)  % loop over number of mitochondrial images
    
    % Setting up for ROI file list
    mname = files(l).name;  % grab mito image name
    mpath = files(l).folder; % grab mito image folder
    flist = get_all_files('roi',[pwd,filesep,mname(1),'-Roi']); % get all roi's for present working directory
    
    for i = 1:numel(flist) % loop over roi files to get them into matlab working memory
        sroi(i) = getimjroi([flist(i).folder,filesep,flist(i).name]);
    end
    
    i1 = fitsread([mpath,filesep,mname]); % read mito image
    i2 = mean(i1,3); % calculate average pixel values
    clear i1 % clear image stack
    
    ils = sort(i2(:));  % sort from lowest to highest pixel
    [m,n] = size(i2); % get the size of image
    
    bkgn = mean(ils(1:round(0.2*numel(ils)))); % background is defined as
    figure
    imagesc(i2);  % show image
    hold on
    A{1,1} = 'ROI Number';
    A{1,2} = 'Mito Ratio';
    for i = 1:numel(sroi)  % loop over all regions of interest
        % Grab pixel region
        top = max(m-sroi(i).vnRectBounds([1,3])+1); % make corrections from imageJ
        left = min(sroi(i).vnRectBounds([2,4])+1);
        bottom = min(m-sroi(i).vnRectBounds([1,3])+1);
        right = max(sroi(i).vnRectBounds([2,4])+1);
        nrat(i) = mean(mean(i2(bottom:top,left:right)))/bkgn; % calculate ration which is average pixel value in region / bkgn
        A{i+1,1} = str2num(sroi(i).strName(4:end)); % prepare data for spreadsheet
        A{i+1,2} =  nrat(i);
        plot([left, right, right, left, left],[bottom, bottom, top, top, bottom],'r'); % plot roi box
    end
axis image
    xlswrite(['cell_',num2str(mname(1)),'_results.xlsx'],A); % write spreadsheet for cell
    totrat = [totrat;nrat.'];
    clear nrat A sroi
    hold off
end

% Display all calcualted ratios
figure
plot(totrat,'.');
xlabel('Index Number');
ylabel('Ratio Value')
title('Mitochondrial Ratios');

% Threshold analysis
for i = 1:100
    thrsh(i) = i/100*(max(totrat)-min(totrat))+min(totrat);
    ind = totrat >= thrsh(i);
    nums(i) = sum(ind);
end
% display threshold results
figure
plot(thrsh,nums)
xlabel('Thershold');
ylabel('Number of regions that pass threshold');
title('Threshold Plot');