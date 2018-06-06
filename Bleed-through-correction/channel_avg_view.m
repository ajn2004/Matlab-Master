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
if ~exist('AVG') % If a correction folder doesn't already exist
    mkdir AVG % make one
end

o = numel(finf);  % number of fits files found

% Bleed Through Correction
for i = 1:o
    fnm = [fpath,finf(i).name]; % concat file name
    imag = fitsinfo(fnm); % get file info
    i1 = fitsread(fnm,'Info', imag); % get file data
    [m,n,p] = size(i1); % get size of data
    
    evens = 2:2:p;
    odds = 1:2:p;
    ie = mean(i1(:,:,evens),3);
    io = mean(i1(:,:,odds),3);
    figure('Units','Normalized','OuterPosition', [0 0 1 1]);
    subplot(1,2,1);
    imagesc(io);
    title('Average of Odd Channels');
    subplot(1,2,2);
    imagesc(ie);
    title('Average of Even Channels');
    
    write2tiff(uint16(ie),['AVG\',finf(i).name(1:end-5),'_avg_even']); % write the red file to the BTC folder to keep separated from raw data
    write2tiff(uint16(io),['AVG\',finf(i).name(1:end-5),'_avg_odd']); % write the green file to the BTC folder to keep separated from raw data
end



%% Local Functions
function write2tiff(i1, fname)
count = 2;
% disp('Yeah it is working');
[m,n,p] = size(i1);
if ~strcmp(fname(end-3:end),'tif')
    fname = [fname,'.tif'];
end
imwrite(uint16(i1(:,:,1)),fname)
if p >1
    while true
        try
            imwrite(uint16(i1(:,:,count)),fname,'WriteMode','append');
            count = count +1;
            if count > numel(i1(1,1,:))
                break
            end
        catch lsterr
            lsterr;
        end
    end
end
end

