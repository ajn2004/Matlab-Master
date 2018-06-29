% Tropical Storm
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
clearvars; close all; clc;

%% USER VARIABLES
pscans = 75;  % Number of prescans as told by labview
fps = 25; % Frames per stimulus: equal to 'Stim Every' value in labview
q = 0.133;  % um/ pixel must be measured for experimental setup
ave_fms = 5; % number of frames pre-stimulus to average to get background
imsafter = 1; % number of images after the stimulus to observe
pixw = 4;  % Window for cutting out region

%% END USER INVOLVEMENT
[fname, fpath] = uigetfile('*local*');  % so far all localization experiments have been named with local in file name
cd(fpath);
c = scrub_config(); % get imaging information.
pix2pho = em_gain(c.Gain);
try % attempt to load dark current info
    load('back_subtract.mat')
catch lsterr
    mi1 = 0;
end
wind = -pixw:pixw; % array variable for selecting square region for fitting

%% Image Loading
i1 = (readtiff(fname) - mi1)/pix2pho; % load image, subtract dark current and convert to photons
% ip1 = rollingball(i1); % rolling ball background subtraction of all frames
ip1 = i1;  % attempt to change up workflow to subtract stimulus response before rolling ball subtraction
[~,~,o] = size(ip1); % grab size of images
TF = sum(sum(i1));
TF =TF(:);
plot(TF);
figure
clear i1
% This next section goes to each frame corresponding to a stimulus,
% averages the previous ave_fms number of frames and subtracts that from
% the subsequent imsafter frames. These background subtracted frames are
% then concatenated into a single variable dip1 and the corresponding
% framenumbers are saved in the fms variable

nstim = floor((o-pscans)/fps); % number of stims
dip1 = [];
fms = [];

for i = 1:nstim % loop over stimuli
    ind = pscans + i*fps; % index now equals stimulus frame
    dip1 = cat(3,dip1, ip1(:,:,ind:ind+imsafter) - mean(ip1(:,:,ind-1-ave_fms:ind-1),3));
    fms =[fms,ind:ind + imsafter];
end
%% remove negative values
dip1 = (dip1 > 0).*dip1;
% dip1 = ip1;
% clear ip1
% M = [];
[~,~,o] = size(dip1);
dip1 = rollingball(dip1); % background subtraction
thrsh = 0.8*mean(max(max(dip1)));
dps = get_das_peaks(dip1,thrsh); % peak detection
[ilocs, fnum, cents] = divide_up(dip1, pixw, dps); % image segmentation
for i = 1:o
    imagesc(dip1(:,:,i))
    ind = fnum == i;
    draw_boxes(cents(ind,:),pixw);
%     title(['Relative Frame ', num2str((i)),' absolute frame ', num2str(fms(i))]);
    axis image
    M(i) = getframe(gcf);
end
imagesc(std(dip1,1,3));
% secs = input('Number of sections to analyze? ');
% ilocs = [];
% cents = [];
% fnum = [];
% for i = 1:secs
%     [x,y] = ginput(1);
%     boxes([x,y],pixw,'r');
%     ilocs = cat(3,ilocs,dip1(round(y) + wind, round(x) + wind,:));
%     cen = [x,y];
%     cents = [cents;repmat(cen,o,1)];
%     fnum = [fnum;fms];
% end

% tic
% thrsh = 3*std(iprod(:)) + mean(iprod(:));
% toc
% thrsh = thresh/100*mean(max(max(diprod)));


[ibkgn] = chop_up(ip1,pixw,fms(fnum), cents,imsafter, ave_fms);
% ilocs = reshape(ilocs,numel(ilocs(:,1,1))^2,numel(ilocs(1,1,:)));
[snr, sig, nois] = sig2noi(ilocs, ibkgn,5);
histogram(sig)
hold on
histogram(nois);
hold off
figure
% % localize selected regions
cal = load('bead_astig_3dcal.mat');
[xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(ilocs, fnum, cents, cal, pixw);
zf_all = zf_all/q;
zf_crlb = zf_crlb./q^2;
Points_diag;
% [xf_all,xf_crlb, yf_all,yf_crlb,sigx_all, sigx_crlb, sigy_all, sigy_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs_sigs(ilocs, fnum, cents, 0);