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
q = 0.128;  % um/ pixel must be measured for experimental setup
pixw = 4;  % Window for cutting out region
fps = 6;    % frames per set is the number of frames / tiff stack
bk_fms = 3;    % number of frames of background for measurement
%% END USER INVOLVEMENT
[fname, fpath] = uigetfile('*local*');  % so far all localization experiments have been named with local in file name
cd(fpath);
% c = scrub_config(); % get imaging information.
pix2pho = 1;
% pix2pho = em_gain(300);
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
% [iprod,ip2] = rollingball(i1); % rolling ball background subtraction of all frames
ip1 = i1;
[m,n,o] = size(ip1); % grab size of images

clear i1
% This next section goes to each frame corresponding to a stimulus,
% averages the previous ave_fms number of frames and subtracts that from
% the subsequent imsafter frames. These background subtracted frames are
% then concatenated into a single variable dip1 and the corresponding
% framenumbers are saved in the fms variable
mins = (0.61*0.525/1.4)/(2*q);
maxs = 2*mins;

% sigma = 1.2;
% iprod = rollingball(ip1);
bkg = zeros(m,n);
for i = 1:bk_fms
    bkg = bkg + imgaussfilt(ip1(:,:,i),mins);
end
bkg = bkg / bk_fms;
dip1(:,:,1) = ip1(:,:,bk_fms+1) - bkg;
fms = bk_fms+1;
for j = 1:fps-bk_fms - 1
    dip1(:,:,j+1) = ip1(:,:,bk_fms + j+1) - bkg;
    fms = [fms, bk_fms + j+1];
end
% dsi1(:,:,1) = movsum(ip1(:,:,2),3) - movsum(ip1(:,:,1),3);
% dsi1(:,:,2) = movsum(ip1(:,:,3),3) - movsum(ip1(:,:,1),3);
% fms = [];
stimdex = [];
for i = 1:(o-fps)/fps % loop over stimuli
    tic
    bkg = bkg*0;
    ind = i*fps+bk_fms+1; % index now equals stimulus frame
    for j = 1:bk_fms
        bkg = bkg + imgaussfilt(ip1(:,:,ind - j),mins);
    end
    bkg = bkg / bk_fms;
    
    dip1 = cat(3,dip1, ip1(:,:,ind) - bkg);   % grab stim1 frame
    fms = [fms,ind]; % keep track of frame molecule was found on
    stimdex = [stimdex,ind];
    for j = 1:fps - bk_fms-1
        dip1 = cat(3,dip1, ip1(:,:,ind+j) - bkg); % grab after stim frames
        fms = [fms,ind + j]; % keep track of frame molecule may appear on
    end
end
rip1 = roball(dip1,6,4); % rolling ball the raw images
% rip1 = dip1;
dip2 = bandpass(dip1,mins,maxs);
% dip2 = denoise_psf(rip1,2); % wavelet decomposition w/ watershed threshold @ 2xstd of the first wavelet plane

% dps = get_das_peaks(wip1,2);
%% remove negative values
[m2,n2,o2] = size(dip1);
surf(max(dip2,[],3));
thrsh = input('What is the threshold?');
dps = cpu_peaks(dip2,thrsh,pixw);
for i = 1:o2
    imagesc(rip1(:,:,i));
    %     drawnow
    colormap('gray');
    hold on
    [row,col] = find(dps(:,:,i) == 1);
    %     plot(col,row,'rx');
    
    draw_boxes([col,row],pixw);
    axis image
    drawnow;
    hold off
    %     mvy(i) = getframe(gcf);
    %     waitforbuttonpress
end

[sdi1, fnum, cents] = divide_up(rip1,pixw, dps);


load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat');
%% Fitting Section of code

% Preallocate variables
xf = [];
yf = [];
N = [];
sx = [];
sy = [];
O = [];
xfc = [];
yfx = [];
Oc =[];
Nc = [];
sxc = [];
syc = [];
yfc = [];
llv = [];
fnumb = [];
fnumb = [];
[fits,crlb,lv, fnout] = slim_locs(sdi1,fnum,cents,cal.ang,50,100); % actual fit
% Assign fits to fitting variables
fnout = fnout.';                                                        % save framenumber relative to difference set (to get absolute frame number use fms(fnout)
xf = [xf;fits(:,1)];                                                    % X-Position
yf = [yf;fits(:,2)];                                                    % Y-Position
N = [N; fits(:,3)];                                                     % Number of Photons
sx = [sx;fits(:,4)];                                                    % Sigma in x' direction
sy = [sy;fits(:,5)];                                                    % Sigma in y' direction
O = [O; fits(:,6)];                                                     % Offset
% Lower bound on Variance of fitted variables
xfc = [xfc;crlb(:,1)];
yfc = [yfc;crlb(:,2)];
Nc = [Nc;crlb(:,3)];
Oc = [Oc;crlb(:,6)];
sxc = [sxc;crlb(:,4)];
syc = [syc;crlb(:,5)];
llv = [llv;-abs(lv)];                                                   % Log Likelihood Value
fnumb = [fnumb;fnout];                                          % Correct the Frame number based off correlation result
zf = getdz(sx,sy,cal.z_cal)/q;  % get z values from sigma measurements
coords = [fits(:,1:2),zf];
[ncoords] = astig_tilt(coords,cal);
%     icoords = [xf,yf,zf]; % icoords will contain the initial 'uncorrected' coordinates
%     clear xf yf zf
%     [coords] = astig_tilt(icoords,cal); % Correct tilt induced issues
%     xf = coords(:,1); % Reassign coordinates based off of corrected values
%     yf = coords(:,2);
%     zf = coords(:,3);

% end
% lv = lv.';
ind = N > 0 & N < 3000;
histogram(lv(ind)./N(ind));
t = input('What Threshold?');

ind = ind & lv./N > t;
ind = ind & abs(sx)*2 > 1 & abs(sx) *2 <10;
ind = ind & sy*2 > 1 & sy *2 < 20;

ind = ind & lv./N > t;
ind = ind & sx*2 > 1.5 & sx *2 < 6;
ind = ind & sy*2 > 1.5 & sy *2 < 6;
ind = ind & zf*q < 0.6 & zf*q > -0.6;
imagesc(mean(rip1,3))
hold on
plot(xf(ind),yf(ind),'.')
[x,y] = ginput(2);
ind = ind & xf < max(x) & xf > min(x);
ind = ind & yf < max(y) & yf > min(y);

ind = ind & lv./N > -1.5;
ind = ind & abs(sx)*2 > 1.5 & abs(sx) *2 <10;
ind = ind & sy*2 > 1.5 & sy *2 < 20;

ind = ind & xfc.^0.5*q < 0.1 & yfc.^0.5*q < 0.1;
s0 = (sy.*sx).^0.5;
lp2 = ((q*s0).^2+q^2/12)./N + 8*pi*(q*s0).^4.*O./(q^2*N.^2);
lp = lp2.^0.5;


% zf_all = zf_all/q;
% zf_crlb = zf_crlb./q^2;
% [~,~,o] = size(sdi1);
% fluor = reshape(sum(sum(sdi1)),o,1);
% ind = N < fluor;
% ind = ind & zf_all*q < 0.5 & zf_all > -0.5;
% ind = ind & q*xf_crlb.^0.5 < 1 & q*yf_crlb.^0.5 < 1;
wave_trajectories;
% figure
xs = xf(ind)*q;
ys = xf(ind)*q;
zs = xf(ind)*q;
% fs = fnum
frate = sum(ind)/numel(ind);
save('Analysis.mat');