%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bracket Quhzx Neural 1
%
% Created AJN 5-6-14
% Updated AJN 11/4/15
% This script will analyze a folder of tiff's and allows you to choose an
% initial threshold and background. The script batch localizes the folder
% chosen, and can be modified to systematically investigate changes in
% variables of interest (such as threshold or rolling ball radius)
% After backgrounds and thresholds have been selected, a batch_levels.m
% file is created saving these values for use later

% This will work for either 1 color or multi-color localization by
% selecting the appropriate channel_number(e.g. 1 for single 2 for multi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc; close all;
p = mfilename('fullpath');
[fpath, fname, fext] = fileparts([p, '.m']);
%% Variables and User modified section
channel_number = 1;                               % Corresponds to the number of channels on your chip, will determine type of localization
an_dir = 'C:\Users\AJN Lab\Desktop\5-30-17 munc13-halo\Subfolder\';   %Correlation file needs to be in this same folder
% an_dir = 0;
corr_name='corr_data.mat';                  % for multi-color must have correlation file in the analysis directory
% email_address = {'andrew.j.nelson@maine.edu','lisa.weatherly@maine.edu'};
email_address = {'andrew.j.nelson@maine.edu'};
% ,'lisa.weatherly@maine.edu'};
wvlngth = 549;                                   %  peak wavelength of fluorophore emission
NA = 1.4;
q = 0.39;             % Pixel Size
n_start = 1;        % Starting frame
n_end = 0;             % End 0 for end frame set to 0 for final frame
n_bkgn = 100;          % Frame to use for background measurement
pix2pho = 33.6;          % Pixel to photon ratio (for sCMOS set to 1 after pre-processing)
file_start = 1;        % File to start on
file_end = 0;          % File to end on (Leave value to 0 to analyze all files after file_start)



%% Initialization
[fname1, data_dir] = uigetfile('*.tif', 'Choose a file in the folder');      % Selects File
cd(data_dir);                                                               % Changes active directory to data directory
File_Info=dir('*.tif');                                                     % Lists all tif files in data directory (will include beam profiles and scale images, ensure that only data is in the folder of interest)
nfiles=length(File_Info);                                                   % Number of files being analyzed
start_bkgns = zeros(numel(File_Info),1);
addpath(fpath);
ninfo = dir('Neural_thetas_it_*');

% load(ninfo(1).name);
if file_end == 0
    file_end = nfiles;
end

% gpud = gpuDevice;
amem = 5.9750*10^09;
if channel_number == 1                  
%% Single Color Batching
%     %% Single Color Setup
%     if exist('batch_levels.mat')==2         % Checks for batch_levels.mat file if already created
%         che = input('A batch level file has been found, would you like to use it? ','s');
%         if strcmp(che,'y')||strcmp(che,'Y')
%             load('batch_levels.mat');
%         else
%             for i=file_start:file_end        
%                 base_name=File_Info(i).name;
%                 [start_bkgns(i)]= func_einzel01_bracket_setup_neural(data_dir,base_name,an_dir,q,n_start,n_bkgn,pix2pho);
%                 
%             end
%             save('batch_levels.mat','start_bkgns','File_Info');
%         end
%     
%     else
%         for i=file_start:file_end
%             base_name=File_Info(i).name;
%             [start_bkgns(i)]= func_einzel01_bracket_setup_neural(data_dir,base_name,an_dir,q,n_start,n_bkgn,pix2pho);
%             
%         end
%         save('batch_levels.mat','start_bkgns','File_Info');
%     end
    
    %% Single Color analysis
%     for j=1:numel(ninfo)  %% andrew modified
        load(ninfo.name);
        gauss_scale = 1;
        bkgn = start_bkgns(1);
        for i = 1:numel(File_Info)
            base_name=File_Info(i).name;
        %           base_name=base_name(1:end-4);
        tic
        func_Quhzx_02_6(data_dir,base_name,an_dir,1,Theta1, Theta2,q,pix2pho,n_start,n_end, wvlngth, NA, n_bkgn, ninfo.name, amem);
        toc
        end
%     end

%% 2 color setup
elseif channel_number == 2
%% Multi-color Batching
    %% Multi color Setup
    if exist('batch_levels.mat')==2
     che = input('A batch level file has been found, would you like to use it? ','s');
     if strcmp(che,'y')||strcmp(che,'Y')
         load('batch_levels.mat');
     else
        for i=1:nfiles        
            base_name=File_Info(i).name;

            [start_bkgns(i)]= func_einzel_2color05_setup(data_dir,base_name,n_bkgn,pix_to_pho);
         end
            save('batch_levels.mat','start_bkgns','File_Info');
     end
    else

         for i=1:nfiles   
              base_name=File_Info(i).name;
              [start_bkgns(i)]= func_einzel_2color05_setup(data_dir,base_name,n_bkgn,pix_to_pho);
         end
         save('batch_levels.mat','start_bkgns','File_Info');
    end
    
    %% Multi-color Analysis
    for j=1:nfiles
        gauss_scale = 1;
        bkgn = start_bkgns(j);
        base_name=File_Info(j).name;
        %           base_name=base_name(1:end-4);

        func_Quhzx_03_1(data_dir,base_name,an_dir,bkgn,Theta1, Theta2,q,pix2pho,n_start,n_end, wvlngth, NA, n_bkgn,corr_name);

    end
end

[reg, name] = system('hostname');
message = ['Your analysis has completed'];
subject = ['Analysis complete'];
matlabmail(email_address, message, subject , 'matlab4ajn@gmail.com','matlabrules');