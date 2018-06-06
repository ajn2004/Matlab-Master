%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organize_Data.m
%
% This script will organize data generated from region_measure_2.m
%
% AJN 9/22/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

nfiles = 1;

% Use a try catch loop to load as many files as necessary, handles variable
% number of files
try
    while true
        
        [fname_temp, fpath_temp] = uigetfile('*mat','Select a file in chronological order, cancel to continue');
        cd(fpath_temp);
        fname{nfiles} = cellstr(fname_temp);
        fpath{nfiles} = cellstr(fpath_temp);
        nfiles = nfiles +1;
    end
catch lasterr;
    nfiles=nfiles - 1;
    disp(['Combining ', num2str(nfiles), ' files']);
end

%% Load structures and add them to total data structure
for i = 1:nfiles
    load([char(fpath{i}), char(fname{i})], 's', 'area_ind', 'grid_size');
    q = 0.128;
    if i == 1
        s_tot = s(area_ind);
        grids = ones(numel(area_ind),1)*grid_size;
        qs = ones(numel(area_ind),1)*q;
    else
        s_tot = vertcat(s_tot, s(area_ind));
        grids = vertcat(grids, ones(numel(area_ind),1)*grid_size);
        qs = vertcat(qs, ones(numel(area_ind), 1)*q);
    end
    clear s area_ind q grid_size
end

%% parse total data structure to perform stats easier
fields = fieldnames(s_tot);
for i = 1:numel(fields); % this loop creates variables for all single vector fields in s
    v = genvarname(strcat(fields(i))); %generate variable name based on relevant field name
    try
        eval(char(strcat(v, ' = cat(1, s_tot.(fields{i}));'))); %forces matlab to evaluate this string
        if strcmp(v{1},'Area')
            eval(char(strcat('men = mean(',v{1},'.*grids.^2./1000000)')));
            eval(char(strcat('stan = std(',v{1},'.*grids.^2./1000000)')));
            disp(['The mean value of the ', v{1},' is ', num2str(men),' +/- ', num2str(stan),' um^2']);
        end
        if strcmp(v{1} ,'Perimeter') || strcmp(v{1}, 'MajorAxisLength') || strcmp(v{1}, 'MinorAxisLength')
            eval(char(strcat('men = mean(',v{1},'.*grids./1000);')));
            eval(char(strcat('stan = std(',v{1},'.*grids./1000);')));
            disp(['The mean value of the ', v{1},' is ', num2str(men),' +/- ', num2str(stan), ' um']);
        end
        if strcmp(v{1}, 'Eccentricity') ||strcmp(v{1}, 'Solidity')
            eval(char(strcat('men = mean(',v{1},');')));
            eval(char(strcat('stan = std(',v{1},');')));
            disp(['The mean value of the ', v{1},' is ', num2str(men),' +/- ', num2str(stan)]);
        end
    catch lasterr
    end
end

%% Represent data
f1 = figure
screen_size = get(0, 'ScreenSize');
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
subplot(5,5,1); plot(Area.*grids.*grids/1000000, MajorAxisLength.*grids./1000, '.b');
xlabel('Area in um^2');
ylabel('MajorAxisLength in um');
subplot(5,5,2); plot(Area.*grids.*grids/1000000, MinorAxisLength.*grids./1000, '.b');
xlabel('Area in um^2');
ylabel('MinorAxisLength in um');
subplot(5,5,3); plot(Area.*grids.*grids/1000000, Perimeter.*grids./1000, '.b');
xlabel('Area in um^2');
ylabel('Perimeter in um');
subplot(5,5,4); plot(Area.*grids.*grids/1000000, Solidity, '.b');
xlabel('Area in um^2');
ylabel('Solidity');
subplot(5,5,5); plot(Area.*grids.*grids/1000000, Eccentricity, '.b');
xlabel('Area in um^2');
ylabel('Eccentricity');

subplot(5,5,6); plot(MajorAxisLength.*grids./1000, Area.*grids.*grids/1000000, '.b');
ylabel('Area in um^2');
xlabel('MajorAxisLength in um');
subplot(5,5,7); plot(MajorAxisLength.*grids./1000, MinorAxisLength.*grids./1000, '.b');
xlabel('MajorAxisLength in um');
ylabel('MinorAxisLength in um');
subplot(5,5,8); plot(MajorAxisLength.*grids./1000, Perimeter.*grids./1000, '.b');
xlabel('MajorAxisLength in um')
ylabel('Perimeter in um');
subplot(5,5,9); plot(MajorAxisLength.*grids./1000, Solidity, '.b');
xlabel('MajorAxisLength in um')
ylabel('Solidity');
subplot(5,5,10); plot(MajorAxisLength.*grids./1000, Eccentricity, '.b');
xlabel('MajorAxisLength in um')
ylabel('Eccentricity');

subplot(5,5,11); plot(MinorAxisLength.*grids./1000, Area.*grids.*grids/1000000, '.b');
xlabel('MinorAxisLength in um');
ylabel('Area in um^2');
subplot(5,5,12); plot(MinorAxisLength.*grids./1000, MajorAxisLength.*grids./1000,  '.b');
ylabel('MajorAxisLength in um');
xlabel('MinorAxisLength in um');
subplot(5,5,13); plot(MinorAxisLength.*grids./1000, Perimeter.*grids./1000, '.b');
xlabel('MinorAxisLength in um');
ylabel('Perimeter in um');
subplot(5,5,14); plot(MinorAxisLength.*grids./1000, Solidity, '.b');
xlabel('MinorAxisLength in um');
ylabel('Solidity');
subplot(5,5,15); plot(MinorAxisLength.*grids./1000, Eccentricity, '.b');
xlabel('MinorAxisLength in um');
ylabel('Eccentricity');

subplot(5,5,16); plot(Perimeter.*grids./1000, Area.*grids.*grids/1000000, '.b');
xlabel('Perimeter in um');
ylabel('Area in um^2');
subplot(5,5,17); plot(Perimeter.*grids./1000, MajorAxisLength.*grids./1000,  '.b');
xlabel('Perimeter in um');
ylabel('MinorAxisLength in um');
subplot(5,5,18); plot(Perimeter.*grids./1000, MinorAxisLength.*grids./1000,  '.b');
ylabel('MinorAxisLength in um');
xlabel('Perimeter in um');
subplot(5,5,19); plot(Perimeter.*grids./1000, Solidity, '.b');
xlabel('Perimeter in um');
ylabel('Solidity');
subplot(5,5,20); plot(Perimeter.*grids./1000, Eccentricity, '.b');
xlabel('Perimeter in um');
ylabel('Eccentricity');

subplot(5,5,21); plot( Solidity, Area.*grids.*grids/1000000, '.b');
xlabel('Solidity');
ylabel('Area in um^2');
subplot(5,5,22); plot( Solidity, MajorAxisLength.*grids./1000,  '.b');
xlabel('Solidity');
ylabel('MinorAxisLength in um');
subplot(5,5,23); plot( Solidity, MinorAxisLength.*grids./1000,  '.b');
xlabel('Solidity');
ylabel('MInorAxisLength in um');
subplot(5,5,24); plot( Solidity, Perimeter.*grids./1000, '.b');
ylabel('Perimeter in um');
xlabel('Solidity');
subplot(5,5,25); plot( Solidity, Eccentricity, '.b');
xlabel('Solidity');
ylabel('Eccentricity');
