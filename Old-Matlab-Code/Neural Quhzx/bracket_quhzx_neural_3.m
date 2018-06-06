%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bracket Quhzx Neural 3
%
% Created AJN 5-6-14
% Updated AJN 11/4/15
% This script will analyze a folder of tiff's and allows you to choose an
% initial threshold and background. The script batch localizes the folder
% chosen, and can be modified to systematically investigate changes in
% variables of interest (such as threshold or rolling ball radius)
% After backgrounds and thresholds have been selected, a batch_levels.m
% file is created saving these values for use later

% This version only supports 1 channel analysis but incorporates fourier
% background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc; close all;
p = mfilename('fullpath');
[fpath, fname, fext] = fileparts([p, '.m']);
%% Variables and User modified section
channel_number = 1;                               % Corresponds to the number of channels on your chip, will determine type of localization
an_dir = 'C:\Users\AJN Lab\Desktop\5-30-17 munc13-halo\Subfolder\';   %Correlation file needs to be in this same folder
corr_name='corr_data.mat';                  % for multi-color must have correlation file in the analysis directory
email_address = {'andrew.j.nelson@maine.edu'};

sp = 50;
wvlngth = 549;                                   %  peak wavelength of fluorophore emission
NA = 1.4;
q = 0.39;             % Pixel Size
n_start = 1;        % Starting frame
n_end = 0;             % End 0 for end frame set to 0 for final frame
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

%% Single Color Batching

%% Single Color analysis
load(ninfo.name);
gauss_scale = 1;
bkgn = start_bkgns(1);
for i = 1:numel(File_Info)
    base_name=File_Info(i).name;
    func_Quhzx_02_7(data_dir,base_name,an_dir,1,Theta1, Theta2,q,pix2pho,n_start,n_end, wvlngth, NA, ninfo.name, amem, sp);
end


[reg, name] = system('hostname');
message = ['Your analysis has completed'];
subject = ['Analysis complete'];
matlabmail(email_address, message, subject , 'matlab4ajn@gmail.com','matlabrules');