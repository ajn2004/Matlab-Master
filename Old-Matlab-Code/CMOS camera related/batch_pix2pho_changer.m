%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to convert sCMOS data so that all pixels have a pixel to photon
% ratio of 1
% This is to be used after pix2photo_conversion.m
% AJN 2/12/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Load pixel to photon map
[fname1 fpath1] = uigetfile('pix2photo_map.mat','Choose the Pixel to Photon Map data set.');
load([fpath1, fname1],'pixmap');
    
%% Choose image folder to analyze
[fname2 fpath2] = uigetfile('*.tif','Choose an image file');
cd(fpath2);
dir_ind= dir('*.tif');

%% For loop to run over all files in folder
for filenum = 1:numel(dir_ind.name)
    clear imagfo pho1 im1 lki
    
    %% If structure to break down file size automatically
    if dir_ind(filenum).bytes <= 4*10^8;                    %if smaller than 400 MB treat the file a whole
        imagfo= imfinfo(dir_ind(filenum).name);
        pho1 = zeros(imagfo(1).Height,imagfo(1).Width);
        for frm_num = 1:numel(imagfo)
        	pho1(:,:,frm_num) = double(imread(dir_ind(filenum).name,frm_num))./pixmap;
        end
        save([dir_ind(filenum).name(1:length(dir_ind(filenum).name)-3), 'mat'],'pho1'); 
    else                                                    % If larger than 400 MB automatically split up the file
        check = 0;  
        count = 1;
        %% While loop automatically finds ideal chunk size count is the number of chunks
        while check == 0       
            dat_size = dir_ind(filenum).bytes/count;
            if dat_size <= 4*10^8
                check = 1;
            else
                count = count+1;
            end
        end
        
        imagfo= imfinfo(dir_ind(filenum).name);             % information about the image
        chunk_im = floor(numel(imagfo)/count);              % chunk size to be analyzed
        
        

        for chunk_num = 1:count
            clear im1 pho1 lki
            strt_frm = (chunk_num-1)*chunk_im+1;
            end_frm  = (chunk_num-1)*chunk_im+chunk_im;
            dum_ind =1;
            pho1 = zeros(imagfo(1).Height,imagfo(1).Width);
            for frm_num = strt_frm:end_frm
                pho1(:,:,dum_ind) = double(imread(dir_ind(filenum).name,frm_num))./pixmap;
                dum_ind = dum_ind + 1;
            end
            save([dir_ind(filenum).name(1:length(dir_ind(filenum).name)-4), '_frames_', num2str(strt_frm), '_', num2str(end_frm), '.mat'],'pho1'); 
        end
         %% Check for a remaineder   
        remain = numel(imagfo)-count*chunk_im;
        if remain ~= 0
           clear pho1
           dum_ind = 1;
           pho1 = zeros(imagfo(1).Height,imagfo(1).Width);
           for framin = end_frm+1:end_frm+1+remain
               pho1(:,:,dum_ind) = double(imread(dir_ind(filenum).name,framin))./pixmap;
               dum_ind = dum_ind+1;
           end
           save([dir_ind(filenum).name(1:length(dir_ind(filenum).name)-4), '_frames_', num2str(end_frm), '_', num2str(end_frm+1+remain), '.mat'],'pho1');
        end
    end
end

clear all