close all; clearvars; clc;  % Standard clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Massive Bleedthrough Correction
%
% A simple script that will allow the user to identify the first or second
% frame as the 'red' channel, then multiply that frame by a user specified
% value and subtract it from the 'green' channel. The 'massive' version
% selects all files in a folder and subfolder
%
% AJN 6/27/18 : Ryan Lab  @ ajn2004@med.cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER VARIABLES
reds = 2;    % Frame for first red channel
bt = 21;      % bleed through in % (i.e. for 20% enter 20)
st_fl = 'C:\Users\AJN Lab\Dropbox\Data\Ryan'; % starting folder to sort through
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finf = get_all_files('fits',st_fl); % grab all fits files in subfolders

o = numel(finf);  % number of fits files found


% Bleed Through Correction
for i = 1:o
    fnm = [finf(i).folder,'\',finf(i).name]; % concat file name
    if ~contains(fnm,'BTC')
        imag = fitsinfo(fnm); % get file info
        i1 = fitsread(fnm); % get file data
        [m,n,p] = size(i1); % get size of data
        
        imaxs = i1(:,:,reds:2:p); % separate out the red channel from the green
        
        if reds == 2 % if the red channel is frame 2
            iblds = i1(:,:,reds - 1:2:p); % green channels are odd
        else % if the red channel is frame 1
            iblds = i1(:,:,reds + 1:2:p); % green channels are even
        end
        
        if reds == 2
            ig = i1(:,:,reds-1:2:p) - (bt/100) * i1(:,:,reds:2:p);
        elseif reds(i) == 1
            ig = i1(:,:,reds+1:2:p) - (bt/100) * i1(:,:,reds:2:p);
        end
        mkdir([finf(i).folder,'\BTC']);
        writetiff(uint16(imaxs),[finf(i).folder,'\BTC\',finf(i).name(1:end-5),'_btc_red']); % write the red file to the BTC folder to keep separated from raw data
        writetiff(uint16(ig),[finf(i).folder,'\BTC\',finf(i).name(1:end-5),'_btc_green']); % write the green file to the BTC folder to keep separated from raw data
    end
end

clear vars
close all
clc