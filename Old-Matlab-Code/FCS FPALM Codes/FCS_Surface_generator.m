clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FCS Surface Generator
% This is a script to generate a surface from FCS excitation rate data
% using Intensity as one axis, exposure time as another, and localization
% uncertainty as the third
%
% AJN 8/12/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User defined variables
int_col = 13; % This is the column from the excel data that the intensity is located
count_col = 11; % This is the column from the excel data that the adjusted counts / molecule/sec are located
bkgn_col = 10; % This is the column from the excel data that the avg background counts are located
max_time = 1;   % Max frame exposure in seconds
time_div = 30; % number of divisions of exposure time (correpsonds to number of data points)

% Thompson Larson Webb related Variables
NA = 1.2; % NA of Lens used
wvlngth = 567; % Wavelength of light imaged in nm
q = 133;  % Pixel Size in nm

%% File selection and loading
[fname, fpath] = uigetfile('*.xls', 'Select data file to analyze');  % forces the user to chose a .xls file
mast_file = xlsread([fpath,fname]);   % loads all xls data

% cherry pick relevant data of interest
intensities = mast_file(:,int_col);
counts_per_mol = mast_file(:,count_col);
bkgns = mast_file(:,bkgn_col);

% Remove NaNs
counts_per_mol(isnan(intensities)) = [];
bkgns(isnan(intensities)) = [];
intensities(isnan(intensities)) = [];
% clear data for proper memory management
clear mast_file 

num_entry = numel(intensities); % find number of elements from which data will be taken

% Construction of grids to be used in calculation
[xgrid, ygrid] = meshgrid(max_time/time_div:max_time/time_div:max_time,counts_per_mol);
[nullgrid, bkgn_grid] = meshgrid(1/time_div:max_time/time_div:max_time,bkgns);
clear nullgrid

xax = max_time/time_div:max_time/time_div:max_time;
% build diffraction information
r0 = 0.61 * wvlngth / NA; % 1 / e^2 radius of PSF
s = r0 /2; % e^-2 radius is 2 sigma, half that is std dev of psf required for TLW

%% TLW equation
lu2 = ((s^2+q^2/12)./(ygrid.*xgrid))+((4*pi^0.5*s^3.*(xgrid.*bkgn_grid).^2)./(q.*(ygrid.*xgrid).^2));
lu = lu2.^0.5;  % loc uncertainty

%% Graphical Representation
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
surf(xax,intensities,lu);  
xlabel('Frame exposure time');
ylabel('Intensity in kW/cm^2');
zlabel('Localization Uncertainty in nm')
ax = gca;
ax.YScale = 'log';
% ax.XScale = 'log';

% %% Resetting Axes
% xticks = get(gca, 'XTick');
% xticks(1) = [];
% yticks = get(gca,'YTick');
% yticks(1) = [];
% set(gca,'XTickLabel',[0,xax(xticks)]);
% set(gca, 'YTickLabel', [0, intensities(yticks).']);