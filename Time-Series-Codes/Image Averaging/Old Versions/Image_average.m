% Image Averager
%
% A simple script to take in several time-averaged trials and create 1
% final averaged copy
%
% Ajn 3/7/17
clear all; close all; clc;

[fname, fpath] = uigetfile('*.tif');

cd(fpath);
files = dir('*.tif');
for j = 1:numel(files)
    finfo = dir([files(j).name(1:end-4),'*']);
    if numel(finfo) > 1
        imagf1 = imfinfo(files(j).name);
        ifin = zeros(imagf1(1).Height,imagf1(1).Width,numel(imagf1));
        for i = 1:numel(finfo)
            imagfo = imfinfo(finfo(i).name);
            for k = 1:numel(imagfo)
                i1(:,:,k) = double(imread(finfo(i).name,k));
            end
            ifin = ifin +i1;
        end
        ifin = round(ifin./numel(finfo));
        imwrite(uint16(ifin(:,:,1)),['Averaged_',fname]);
        for i = 2:numel(imagf1)
            imwrite(uint16(ifin(:,:,i)),['Averaged_',files(j).name],'writemode','append');
        end
        
        sums = sum(sum(ifin,1),2);
        peak = find(sums == max(sums));
        
        avefin = mean(ifin,3);
        peakfin = ifin(:,:,peak) - avefin;
        imwrite(uint16(peakfin), ['Sub_',files(j).name]);
    end
end