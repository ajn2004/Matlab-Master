close all; clearvars; clc;  % Standard clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jaime's Bleedthrough Correction
%
% A simple script that will allow the user to identify the first or second
% frame as the 'red' channel, then multiply that frame by a user specified
% value and subtract it from the 'green' channel
%
% AJN 8/1/17 : Ryan Lab  @ ajn2004@med.cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% File Selection and Starting up
[fname, fpath] = uigetfile('*.fits','Choose a fits file'); % Select a fits file
cd(fpath);  % change matlab's path to the data folder

finf = dir('*.fits'); % find all fits images in the folder
if ~exist('BTC') % If a correction folder doesn't already exist
    mkdir BTC % make one
end

o = numel(finf);  % number of fits files found

%% Building data set for Red Channel Identification
for i = 1:o
    fnm = [fpath,finf(i).name];  % concatenate file name
    imag = fitsinfo(fnm); % get image about file
    i1 = fitsread(fnm,'Info', imag); % load all data into variable i1
    i1t(:,:,:,i) = i1(:,:,1:2); % only select the first 2 frames into a 4D variable
    clear i1 % cleanup
end
%% Set the stage for User ID of Red Channel
close all % close all previously opened figures
figure('Units','Normalized','Outerposition',[0, 0, 1, 1]); % full screen figure
reds  = zeros(o,1); % preallocate array
for i = 1:o % loop over all files
    subplot(1,2,1);imagesc(i1t(:,:,1,i));title('Frame 1');axis image;  % show frame 1
    subplot(1,2,2);imagesc(i1t(:,:,2,i));title('Frame 2');axis image;  % show frame 2
    while true % a while loop to ensure you only enter what I want you to enter (i.e. Failure Analysts eat your hearts out)
        try % Try to get a reasonable answer from the user
            ams = input('Which frame is the red channel?(1/2)'); % ask user for bleedthrough percentage as a number between 0 and 100
        catch lsterr % if user enters a string invoke the following
            clc; % clean up console to not flood with obnoxious lines
            disp('That was not either 1 or 2'); % Chastise the user lightly as they may not know the difference between numbers and letters
        end
        if ams == 1  % if the user entered 1
            reds(i) = 1; % save the value as 1
            break % get out of the loop
        elseif ams == 2 % if the user entered 2
            reds(i) = 2; % save the value as 2
            break % get out of the loop
        else % User entered a different number?
            clc
            disp('Jaime this is not programmed to handle wrong answers, try again'); % Get mad at Jaime!
        end
    end
end

%% almost done with User Input
close all % clean up figures
clear i1t ams i fnm imag % just a little clean up before the next part

% Ask user for bleedthrough percentage as a number between 0 and 100
while true % Ensure the User enters what we want
    try % Try loop incase of a string entered
        bt = input('What percentage bleedthrough did you measure? (11.5 percent would be entered as 11.5)'); % Ask politely and ye shall receive
        break % if that worked hop out
    catch lsterr % if a string was entered
        clc % clean up console
        disp('I need a number'); % Demand a number, don't back down!
    end
end

% Bleed Through Correction
for i = 1:o
    fnm = [fpath,finf(i).name]; % concat file name
    imag = fitsinfo(fnm); % get file info
    i1 = fitsread(fnm,'Info', imag); % get file data
    [m,n,p] = size(i1); % get size of data
    
    imaxs = i1(:,:,reds(i):2:p); % separate out the red channel from the green
    
    if reds(i) == 2 % if the red channel is frame 2
        iblds = i1(:,:,reds(i) - 1:2:p); % green channels are odd
    else % if the red channel is frame 1
        iblds = i1(:,:,reds(i) + 1:2:p); % green channels are even
    end
    
    if reds(i) == 2
        ig = i1(:,:,reds(i)-1:2:p) - (bt/100) * i1(:,:,reds(i):2:p);
    elseif reds(i) == 1
        ig = i1(:,:,reds(i)+1:2:p) - (bt/100) * i1(:,:,reds(i):2:p);
    end
    
    writetiff(uint16(imaxs),['BTC\',finf(i).name(1:end-5),'_btc_red']); % write the red file to the BTC folder to keep separated from raw data
    writetiff(uint16(ig),['BTC\',finf(i).name(1:end-5),'_btc_green']); % write the green file to the BTC folder to keep separated from raw data
end

    % end
    clear ifin i1 cr cg j i icorr iblds imaxs bt
    close all
    clc
    
    
    
    