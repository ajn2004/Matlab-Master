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
molish = 4; % expected number of molecules
q = 0.133;  % um/ pixel must be measured for experimental setup
pixw = 4;  % Window for cutting out region
fps = 3;    % frames per set is the number of frames / tiff stack
%% END USER INVOLVEMENT
[fname, fpath] = uigetfile('*local*');  % so far all localization experiments have been named with local in file name
cd(fpath);
% c = scrub_config(); % get imaging information.
% pix2pho = em_gain(c.Gain);
pix2pho = em_gain(300);
try % attempt to load dark current info
    load('back_subtract.mat')
catch lsterr
    mi1 = 0;
end
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
% ip1 = rollingball(i1); % rolling ball background subtraction of all frames
ip1 = i1;  % attempt to change up workflow to subtract stimulus response before rolling ball subtraction
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
dps = zeros(m,n,2*(o)/3);
    
for i = 1:(o-fps)/3 % loop over stimuli
tic
    ind = i*3+2; % index now equals stimulus frame
    dip1 = cat(3,dip1, ip1(:,:,ind) - ip1(:,:,ind-1));   % grab stim1 frame
    dip1 = cat(3,dip1, ip1(:,:,ind+1) - ip1(:,:,ind-1)); % grab stim2 frame
    fms =[fms,ind,ind+1];
    
    [cents] = eyeindsky(ip1(:,:,ind-1:ind),molish);
    for j = 1:molish
    dps(cents(j,2),cents(j,1),numel(dip1(1,1,:))-1) = 1;
    end
    [cents] = eyeindsky(ip1(:,:,[ind-1,ind+1]),molish);
    for j = 1:molish
    dps(cents(j,2),cents(j,1),numel(dip1(1,1,:))) = 1;
    end
    t(i) = toc;
    ajn_wait(t, i, (o-fps)/3);
end
ip1 = rollingball(ip1); % rolling ball the raw images
%% remove negative values
dip1 = (dip1 > 0).*dip1;

[~,~,o] = size(dip1);
dip1 = rollingball(dip1); % background subtraction




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
%     boxes([x,y],pixw,'r');close all

%     ilocs = cat(3,ilocs,dip1(round(y) + wind, round(x) + wind,:));
%     cen = [x,y];
%     cents = [cents;repmat(cen,o,1)];
%     fnum = [fnum;fms];
% end

% tic
% thrsh = 3*std(iprod(:)) + mean(iprod(:));
% toc
% thrsh = thresh/100*mean(max(max(diprod)));


% [ibkgn] = chop_up(ip1,pixw,fms(fnum), cents,imsafter, ave_fms);
% % ilocs = reshape(ilocs,numel(ilocs(:,1,1))^2,numel(ilocs(1,1,:)));
% [snr, sig, nois] = sig2noi(ilocs, ibkgn,5);
% histogram(sig)
% hold on
% histogram(nois);
% hold off
% figure
% % localize selected regions
cal = load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Single-release-codes\bead_astig_3dcal.mat');
[xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(ilocs, fnum, cents, cal, pixw);
zf_all = zf_all/q;
zf_crlb = zf_crlb./q^2;
% Points_diag;
% [drifts = get_drift_ims(ip1);
Traj_show;
save('Localization_file.mat','xf_all','xf_crlb', 'yf_all','yf_crlb','zf_all', 'zf_crlb', 'N', 'N_crlb','off_all', 'off_crlb', 'framenum_all', 'llv', 'iters');
% [xf_all,xf_crlb, yf_all,yf_crlb,sigx_all, sigx_crlb, sigy_all, sigy_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs_sigs(ilocs, fnum, cents, 0);