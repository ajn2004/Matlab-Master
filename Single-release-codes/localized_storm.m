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
molish = 1; % expected number of molecules
q = 0.133;  % um/ pixel must be measured for experimental setup
pixw = 7;  % Window for cutting out region
fps = 3;    % frames per set is the number of frames / tiff stack
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
mins = (0.61*0.51/1.3)/(2*q);
maxs = 2*mins;

% sigma = 1.2;
% iprod = rollingball(ip1);
dip1(:,:,1) = ip1(:,:,2) - imgaussfilt(ip1(:,:,1),mins);
dip1(:,:,2) = ip1(:,:,3) - imgaussfilt(ip1(:,:,1),mins);
% dsi1(:,:,1) = movsum(ip1(:,:,2),3) - movsum(ip1(:,:,1),3);
% dsi1(:,:,2) = movsum(ip1(:,:,3),3) - movsum(ip1(:,:,1),3);
fms = [];
dps = zeros(m,n,2*(o)/3);
dip2(:,:,1) =     ip1(:,:,2);
dip2(:,:,2) =     ip1(:,:,3);
for i = 1:(o-fps)/3 % loop over stimuli
tic
    ind = i*3+2; % index now equals stimulus frame
    dip1 = cat(3,dip1, ip1(:,:,ind) - imgaussfilt(ip1(:,:,ind-1),mins));   % grab stim1 frame
    dip1 = cat(3,dip1, ip1(:,:,ind+1) - imgaussfilt(ip1(:,:,ind-1),mins)); % grab stim2 frame
%     dip1 = cat(3,dip1, imgaussfilt(ip1(:,:,ind),mins) - imgaussfilt(ip1(:,:,ind-1),mins));   % grab stim1 frame
%     dip1 = cat(3,dip1, imgaussfilt(ip1(:,:,ind+1),mins) - imgaussfilt(ip1(:,:,ind-1),mins)); % grab stim2 frame
%     dsi1 = cat(3,dsi1, movsum(ip1(:,:,ind),3) - movsum(ip1(:,:,ind-1),3));
% 	dsi1 = cat(3,dsi1, movsum(ip1(:,:,ind+1),3) - movsum(ip1(:,:,ind-1),3));
    dip2 = cat(3,dip2, ip1(:,:,ind));
    dip2 = cat(3,dip2, ip1(:,:,ind+1));
    fms =[fms,ind,ind+1];
%     frames  = 
%     [cents] = eyeindsky(ip1(:,:,ind-1:ind),molish);
%     for j = 1:molish
%     dps(cents(j,2),cents(j,1),numel(dip1(1,1,:))-1) = 1;
%     end
%     [cents] = eyeindsky(ip1(:,:,[ind-1,ind+1]),molish);
%     for j = 1:molish
%     dps(cents(j,2),cents(j,1),numel(dip1(1,1,:))) = 1;
%     end
%     t(i) = toc;
%     ajn_wait(t, i, (o-fps)/3);
end
rip1 = rollingball(dip1); % rolling ball the raw images

% dip3 = bandpass(dip1,mins,maxs);
dip2 = denoise_psf(dip1,2); % wavelet decomposition w/ watershed threshold @ 2xstd of the first wavelet plane

% dps = get_das_peaks(wip1,2);
%% remove negative values
[m2,n2,o2] = size(dip1);
% for i = 1:o2
%     dip1(:,:,i) = dip1(:,:,i) - min(min(dip1(:,:,i)));
% end
% dip1 = (dip1 > 0).*dip1;
% dip2 = lp_filt(dip1,4);
% dip2 = bandpass(dip1,0.8,5.5);
surf(max(dip2,[],3));
thrsh = input('What is the threshold?');
% thrsh = max(max(max(dip2(:,:,1:2))));
dps = get_das_peaks(dip2,thrsh);

for i = 1:o2
%     dig(:,:,i) = imgaussfilt(rip1(:,:,i),1.5);
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
% [ind] = find_fm_dupes(cents,fnum,1.5*pixw);
% sdi1(:,:,ind) = [];
% cents(ind,:) = [];
% fnum(ind) = [];

%% Comment section to hold code
% imgaussfilt(dip1,1.5)
%

% [~,~,o] = size(dip1);
% imagesc(sum(dip1,3))
% title('Select an ROI');
% while true 
% [x,y] = ginput(1);
% x = round(x);
% y = round(y);
% draw_boxes([x,y],pixw);
% b = waitforbuttonpress;
% if b == 1
%     break
% end
% hold off
% imagesc(sum(dip1,3));
% end
% wind = -pixw:pixw;
% sdi1 = dip1(y+wind,x+wind,:);
% si1 = ip1(y+wind,x+wind,:);
% rat = [];
% for i = 1:(o-fps)/3
%     ind = i*3+1;
%     rat(numel(rat)+1) = rat_view(si1(:,:,ind:ind+1));
%     rat(numel(rat)+1) = rat_view(si1(:,:,[ind,ind+2]));
% end
% fnum = 1:numel(sdi1(1,1,:));
% <<<<<<< HEAD
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Single-release-codes\z_calib.mat');
% =======
% cal = load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Single-release-codes\bead_astig_3dcal.mat');
% cents = zeros(numel(sdi1(1,1,:)),2);
% [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(sdi1, fnum, cents, cal, pixw);
% >>>>>>> 730d9b17f1dcce149517eda04abbd30f2b8f5031
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
% ang = 0.1047;
% for i = 1:numel(fnum)
   [fits,crlb,lv, fnout] = slim_locs(sdi1,fnum,cents,cal.ang,50,100); 
   fnumb = [];
     fnout = fnout.';                                                        % Save frame number which corresponds to Z-position
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
%     icoords = [xf,yf,zf]; % icoords will contain the initial 'uncorrected' coordinates
%     clear xf yf zf
%     [coords] = astig_tilt(icoords,cal); % Correct tilt induced issues
%     xf = coords(:,1); % Reassign coordinates based off of corrected values
%     yf = coords(:,2);
%     zf = coords(:,3);
    
% end
% lv = lv.';
ind = N > 0 & N < 1500;
histogram(lv(ind)./N(ind));
t = input('What Threshold?');
% <<<<<<< HEAD
ind = ind & lv./N > t;
ind = ind & abs(sx)*2 > 1 & abs(sx) *2 <10;
ind = ind & sy*2 > 1 & sy *2 < 20;
% =======
% <<<<<<< HEAD
ind = ind & lv./N > t;
ind = ind & sx*2 > 1.5 & sx *2 < 6;
ind = ind & sy*2 > 1.5 & sy *2 < 6;
imagesc(mean(rip1,3))
hold on
plot(xf(ind),yf(ind),'.')
[x,y] = ginput(2);
ind = ind & xf < max(x) & xf > min(x);
ind = ind & yf < max(y) & yf > min(y);
% =======
ind = ind & lv./N > -1.5;
ind = ind & abs(sx)*2 > 1.5 & abs(sx) *2 <10;
ind = ind & sy*2 > 1.5 & sy *2 < 20;
% >>>>>>> 730d9b17f1dcce149517eda04abbd30f2b8f5031
ind = ind & xfc.^0.5*q < 0.1 & yfc.^0.5*q < 0.1;
s0 = (sy.*sx).^0.5;
lp2 = ((q*s0).^2+q^2/12)./N + 8*pi*(q*s0).^4.*O./(q^2*N.^2);
lp = lp2.^0.5;

% >>>>>>> 068cdbeba12f439e2ade83caab48adcea5cdf5f6
% zf_all = zf_all/q;
% zf_crlb = zf_crlb./q^2;
% [~,~,o] = size(sdi1);
% fluor = reshape(sum(sum(sdi1)),o,1);
% ind = N < fluor;
% ind = ind & zf_all*q < 0.5 & zf_all > -0.5;
% ind = ind & q*xf_crlb.^0.5 < 1 & q*yf_crlb.^0.5 < 1;
wave_trajectories;
% figure
% plot(fluor);
% xlabel('Framenumber')
% ylabel('Photons')
% title('Sum of every frame');
% figure
% scatter3(xf_all*q, yf_all*q, zf_all*q,[],mod(framenum_all,2));
% axis equal
% zlim([-0.700, 0.4]); % limit of Z in um
% xlabel('Position um')
% ylabel('Position um');
% zlabel('A. Position um')
% title('"Localizations"')
% figure
% % imagesc(std(si1,1,3))
% set_scale(std(sdi1,1,3),q,4);
% colormap('jet')
% hold on
% plot(xf_all(ind),yf_all(ind),'k.')
% axis image
% title('Projection of localizations onto image');
% save('Localization_file.mat','xf_all','xf_crlb', 'yf_all','yf_crlb','zf_all', 'zf_crlb', 'N', 'N_crlb','off_all', 'off_crlb', 'framenum_all', 'llv', 'iters');
% Traj_show;
frate = sum(ind)/numel(ind);
save('Analysis.mat');