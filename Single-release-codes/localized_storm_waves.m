% Localized storm waves
%
% This is a more narrow version of Hurricane designed to analyze data
% associated with the single vesicle release experiments being performed
% currently. The general work flow is as follows, import tiff; subtract
% dark current and convert to photons from ADUs; perform rolling ball
% subtraction step; present user with image representation and ask for cut
% out; send parts of the cut out region to be analyzed based on temporal
% proximity to anticipated stimulus. this will require a configuration file
% be saved in the folder being analyzed
% AJN 6/12/18 Ryan Lab
% localized storm_waves is a new version relying on wavelet analysis
clearvars; close all; clc;

%% USER VARIABLES
q = 0.133;  % um/ pixel must be measured for experimental setup
pixw = 7;  % Window for cutting out region
fps = 3;    % frames per set is the number of frames / tiff stack
%% END USER INVOLVEMENT
% [fname, fpath] = uigetfile('*local*');  % so far all localization experiments have been named with local in file name
% cd(fpath);
% c = scrub_config(); % get imaging information.
% pix2pho = em_gain(c.Gain);
pix2pho = em_gain(300);

wind = -pixw:pixw; % array variable for selecting square region for fitting

%% Image Loading
i1 = [];
files = dir('*tif');
for i = 1:numel(files)
    try
    ind(i) = str2num(files(i).name(7:end-4));
    catch
        ind(i) = 0;
    end
end
[B,I] = sort(ind);

for i = 1:numel(files)
i1 = cat(3,i1,(readtiff(files(I(i)).name))/pix2pho); % load image, subtract dark current and convert to photons
end
% [B,I] = sort(ind); % determine numerical order of files
% [iprod,ip2] = rollingball(i1); % rolling ball background subtraction of all frames
ip1 = i1;
[m,n,o] = size(ip1); % grab size of images

clear i1
% This next section goes to each frame corresponding to a stimulus,
% averages the previous ave_fms number of frames and subtracts that from
% the subsequent imsafter frames. These background subtracted frames are
% then concatenated into a single variable dip1 and the corresponding
% framenumbers are saved in the fms variable
dip1(:,:,1) = ip1(:,:,2) - ip1(:,:,1);
dip1(:,:,2) = ip1(:,:,3) - ip1(:,:,1);
fms = [];
% dps = zeros(m,n,2*(o)/3);
    
for i = 1:(o-fps)/3 % loop over stimuli
tic
    ind = i*3+2; % index now equals stimulus frame
    dip1 = cat(3,dip1, ip1(:,:,ind) - (ip1(:,:,ind-1)));   % grab stim1 frame
    dip1 = cat(3,dip1, ip1(:,:,ind+1) - (ip1(:,:,ind-1))); % grab stim2 frame
    fms =[fms,ind,ind+1];
end

[m2,n2,o2] = size(dip1);
wip1 = denoise_psf(dip1,2); % wavelet decomposition w/ watershed threshold @ 2xstd of the first wavelet plane

dps = get_das_peaks(wip1,2);
[iloc, fnum, cents] = divide_up(wip1,pixw,dps);
wind = -pixw:pixw;
for i = 1:numel(cents(:,1))
    N(i) = sum(sum(dip1(cents(i,2)+wind,cents(i,1)+ wind,fnum(i)).*(dip1(cents(i,2)+wind,cents(i,1)+ wind,fnum(i))>0)));
end
[fits] = wave_fit(iloc);
histogram(N);
xlabel('Sum of Imaged Location');
ylabel('Frequency');
figure
wave_trajectories;