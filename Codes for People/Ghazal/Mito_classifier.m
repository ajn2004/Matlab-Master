clearvars; close all; clc; % clean up
totrat = [];
% [mname, mpath] = uigetfile('*mtRFP*');
files = dir('*mtRFP*');
for l = 1:numel(files)
    mname = files(l).name;
    mpath = files(l).folder;
flist = get_all_files('roi',[pwd,filesep,mname(1),'-Roi']); % get all roi's for present working directory

for i = 1:numel(flist)
    sroi(i) = getimjroi([flist(i).folder,filesep,flist(i).name]);
end


i1 = fitsread([mpath,filesep,mname]);
i2 = mean(i1,3);
clear i1

ils = sort(i2(:));
[m,n] = size(i2);

bkgn = mean(ils(1:round(0.2*numel(ils))));
imagesc(i2);
hold on
A{1,1} = 'ROI Number';
A{1,2} = 'Mito Ratio';
for i = 1:numel(sroi)
    top = max(m-sroi(i).vnRectBounds([1,3]));
    left = min(sroi(i).vnRectBounds([2,4]));
    bottom = min(m -sroi(i).vnRectBounds([1,3]));
    right = max(sroi(i).vnRectBounds([2,4]));
    nrat(i) = sum(sum(i2(bottom:top,left:right)))/bkgn;
    A{i+1,1} = str2num(sroi(i).strName(4:end));
    A{i+1,2} =  nrat(i);
    plot([left, right, right, left, left],[bottom, bottom, top, top, bottom],'r');
end
% plot(nrat)
% 
% xlabel('ROI Number');
% ylabel('Mito over background');
xlswrite(['cell_',num2str(mname(1)),'_results.xlsx'],A);
totrat = [totrat;nrat.'];
end