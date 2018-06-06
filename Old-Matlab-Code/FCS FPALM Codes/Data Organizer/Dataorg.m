%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data organizer
% This is a script that will organize data by intensity and show the
% average values associated with those intensities. It is currently being
% used with the FCS-FPALM Photophysics project
% AJN 10/19/15
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

% Select a file in the directory of interest
[fname, fpath] = uigetfile('*tol.mat','Select a file in folder of interest');
addpath(pwd); % add present working directory to path
cd(fpath); % change working directory to chosen directory
finfo = dir('*tol.mat'); % create a variable of file info
N_tol_max = 250;
lasterror('reset')
try
    load('Organization_params.mat')
    intens = numel(intensities);
catch lasterr
end

if ~strcmp(lasterr,'')
    %% Have user build intensities that were used during the experiment
    intens = input('How many Intensities were tested ');
    intensities = zeros(intens,1);
    for i = 1:intens
        intensities(i) = input(['Value for intensity ', num2str(i), ' ']);
        day(i) = input(['Was this on day 1 or 2? ']);
    end
    
    for i = 1:intens
        range(i,1) = input(['What was the cell number for the first cell at ', num2str(intensities(i)), 'kw/cm^2 on day ', num2str(day(i)), ' ']);
        range(i,2) = input(['What was the cell number for the second cell at ', num2str(intensities(i)), 'kw/cm^2 on day ', num2str(day(i)), ' ']);
    end
    
    save('Organization_params.mat','range','day','intensities')
end
%% Matlab searches the resulting structure finfo to find the starting and ending cells

for i = 1:intens  % loop over number of intensities
    for j = range(i,1):range(i,2) % for each intensity loop over range of cells
        if day(i) == 2
            idx = ~cellfun('isempty',strfind({finfo.name},['cell ', num2str(j)]));  % find all files corresponding to the cell
        else
            idx = ~cellfun('isempty',strfind({finfo.name},['cell', num2str(j)]));  % find all files corresponding to the cell
        end
        index = find(idx == 1); % build index of found files
        for k = 1:numel(index) % loop over files corresponding to the proper index
            clear N bkgn lp
            load(finfo(index(k)).name, 'N','bkgn', 'lp') % load N background noise and LP
            if j == range(i,1) && k == 1  % build vectors over variables of interest
                N_tot = N(N<N_tol_max);
                bkgn_tot = bkgn;
                lp_tot = lp(N<N_tol_max);
            else
                N_tot = vertcat(N_tot,N(N<N_tol_max));
                bkgn_tot = vertcat(bkgn_tot,bkgn);
                lp_tot = vertcat(lp_tot,lp(N<N_tol_max));
            end% finish if statement for building vectors
        end % finish looping over files corresponding to a single cell
    end % finishing looping over cells in a single intensity
    fin_intes(i) = intensities(i);
    Ave_N(i) = mean(N_tot);
    Std_N(i) = std(N_tot);
    Ave_bkgn(i) = mean(bkgn_tot);
    Std_bkgn(i) = std(bkgn_tot);
    Ave_lp(i) = mean(lp_tot);
    Std_lp(i) = std(lp_tot);
end