function [c] = scrub_config()
% Andor's config file is saved as a .cfg format
% Here we'll start a place to reliable get the desired parts of the config
% file and export them as a structured variable c
files = dir('*.cfg');

x = fopen(files(1).name);
y = fread(x);
clear x;

[ind] = find(y == 13); % find every new line return

% Ind now has all a list of every new line where line 1 corresponds to
% y(1:ind(1) and line 120 corresponds to y(ind(119)+1:ind(120)-1)

% get exposure time
sy = y(ind(190)+1:ind(191)-1);
ind1 = find(sy == 61); % Find where there is an equal sign
c.ExposureTime = str2num(char(sy(ind1+1:end)).');
clear ind1 sy

% get Accumulate Cycle Time
sy = y(ind(192)+1:ind(193)-1);
ind1 = find(sy == 61); % Find where there is an equal sign
c.AccumulateCycleTime = str2num(char(sy(ind1+1:end)).');
clear ind1 sy

% get emGain
sy = y(ind(217)+1:ind(218)-1);
ind1 = find(sy == 61); % Find where there is an equal sign
c.Gain = str2num(char(sy(ind1+1:end)).');
clear ind1 sy
