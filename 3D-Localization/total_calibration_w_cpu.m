% Calibration for 3D
% This script is intended to be a one-stop-shop for astigmatic calibration
% of scan data. Comments are heavy to help guide user through programming
% AJN 10-8-18 Ryan Lab

% Preclean UI
clearvars;
close all;
clc;
%% User variables
pixw = 4;                                                                   % ROI 'radius'
q = 0.133;                                                                  % Pixel size in um
step = 10;                                                                  % Steps between frames in nm
CCs = 50;                                                                   % number of frames to x-correlate over
wind = -pixw:pixw;                                                          % create window for segmentation

%% END USER INPUT
% Create display figure as a tab group
f = figure('units','Normalized','OuterPosition',[0 0 1 1]);                 % Initialize figure
tg = uitabgroup(f);                                                         % tg is the tabgroup to help reduce clutter of figures

%% File loading and image segmentation
files = dir('scan*');
t1 = uitab(tg,'Title','Mean Images');                                       % Create tab for image representation
tg1 = uitabgroup(t1);                                                       % Make tabgroup object
A= {};                                                                      % Saving images as a double cell variable
psfs = {};                                                                  % Saving scans of PSF as double cell variable
pind = 1;                                                                   % Counting variable for populating psfs
for i = 1:numel(files)                                                      % Loop over all files
    A{i} = readtiff(files(i).name);                                         % store total images into variable
    ax = axes(uitab(tg1,'Title',files(i).name(1:end-4)));                   % get axes for appropriate tab
    imagesc(ax,max(A{i},[],3))                                              % Represent maximal image
    psf = denoise_psf(max(A{i},[],3),2);                                    % use wavelet transform to identify molecules
    dps = das_peaks(psf,10);                                                % Peak finder
    [row,col] = find(dps == 1);                                             % Find peaks in dps
    ind = find_dupes([col,row],1.5*pixw);                                   % remove overlapping molecules
    % Remove duplicate entries
    row(ind) = [];
    col(ind) = [];
    for j = 1:numel(row)                                                    % loop over all locations, but attempt to center the maximum pixel on each sub image
        draw_boxes([col,row],pixw);                                         % Show location of found regions on representation
        try
            for k = 1:numel(A{i}(1,1,:))
                [r,c] = find(A{i}(row(j) + wind, col(j) + wind,k) == max(max(A{i}(row(j) + wind, col(j) + wind,k)))); % find maxima pixel in subregion   
                psfs{pind}(:,:,k) = A{i}((row(j) - pixw - 1) + r(1) + wind, col(j) - pixw - 1 + c(1) + wind,k); % Store psfs in their own variable
            end
            pind = pind +1;
        catch lsterr
        end
    end
end

% Memory Management
% Here all interested regions are saved in psfs
clear A dps psf row col files pind

%% Cross Correlate PSFs
[m,n,o] = size(psfs{1});                                                    % Grab size of PSF regions
t2 = uitab(tg,'Title','X-Correlation');                                     % Make tab for X-Corr representation
ax = axes(t2);                                                              % Get axes for X-Corr plot

% Correlate all PSFS to the first frame scan
for i = 1:numel(psfs)
    [ind(i),mc] = cross3d(psfs{1},psfs{i}(:,:,round(o/2) + (-CCs:CCs)));    % Perform Correlation
    plot(ax,(1:numel(mc(:))),mc(:)/max(mc(:)))                              % Represent Result
    hold on    
end
hold off
disp = ind - ind(1); % Make Correction to x-corrs by subtracting the identity index 

ang = [];
%% Fitting Analysis
% Get Eliptical Angle
for i = round(o/2) + (-CCs:-10)                                             % Perform eliptical Angle determination over several frames
    [a] = get_elip_ang(psfs{1}(:,:,i),2.5,1.5);                             % Find optimal eliptical angle
    ang = [ang;a];                                                          % Save result to an array
end
ang = mean(ang);                                                            % take mean angle

% Preallocate fitting variables
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
fnum = [];

% CPU fit of PSFs
% for i = 1:numel(psfs)
%     for j = 1:o
%         [y,x] = find(psfs{i}(:,:,j) == max(max(psfs{i}(:,:,j))));
%         [fits,crlb,llv] = func_mle_crlb(psfs{i}(:,:,j),x(1)-pixw-1,y(1)-pixw-1,2,ang);
%     x = fits(1)*cos(-ang) - sin(-ang)*fits(2);
%     y = fits(1)*sin(-ang) + cos(-ang)*fits(2);
%     xf = [xf];
%     yf = [yf];
%     N = [N;fits(3)];
%     O = [O;fits(6)];
%     sx = [sx;abs(fits(4))];
%     sy = [sy;abs(fits(5))];
%     xfc = [xfc;crlb(1)];
%     yfc = [yfc;crlb(2)];
%     Nc = [Nc;crlb(3)];
%     Oc = [Oc;crlb(6)];
%     sxc = [sxc;crlb(4)];
%     syc = [syc;crlb(5)];
%     lv = [lv;-abs(llv)];
%     fnum = [fnum;j - ind(i)];
%     end
% end
% gpu fit of PSFs
fms = 1:o; % input framenum
for i = 1:numel(psfs)
    [fits, crlb, lv,fnout] = slim_locs(psfs{i},fms,zeros(o,2),ang,50,100);
    fnout = fnout.';
    xf = [xf;fits(:,1)];
    yf = [yf;fits(:,2)];
    N = [N; fits(:,3)];
    sx = [sx;fits(:,4)];
    sy = [sy;fits(:,5)];
    O = [O; fits(:,6)];
    xfc = [xfc;crlb(:,1)];
    yfc = [yfc;crlb(:,2)];
    Nc = [Nc;crlb(:,3)];
    Oc = [Oc;crlb(:,6)];
    sxc = [sxc;crlb(:,4)];
    syc = [syc;crlb(:,5)];
    llv = [llv;-abs(lv)];
    fnum = [fnum;fnout - disp(i)]; % given the way the correlation is done, subtracting the correction is the way to go
%     fnum = [fnum;fnout];

end
       z0 = fnum*step; 
       % Tolerance Step
indy = llv./N > -0.1 & N >0 & N < 10000 & abs(syc) < 0.005 & abs(sxc) < 0.005 ;

% Data Representation

t4 = uitab(tg,'Title','Fitting Outputs');
tg4 = uitabgroup(t4);
txy = uitab(tg4,'Title','X-Y Position');
ax = axes(txy);
plot(ax,xf(indy),yf(indy),'.');
axis equal
xlabel('X-Fit')
ylabel('Y-Fit');

tlp = uitab(tg4,'Title','X-Y Uncertainty');
ax = axes(tlp);
histogram(ax,xfc(indy).^0.5*q,'Normalization','Probability');
hold on
histogram(ax,yfc(indy).^0.5*q,'Normalization','Probability');
legend('X-unc','Y-Unc');
xlabel('Uncertainty in um')
ylabel('Probability');
hold off

tsxy = uitab(tg4,'Title','Sigma Uncertainties');
ax = axes(tsxy);
histogram(ax,sxc(indy).^0.5*q,'Normalization','Probability');
hold on
histogram(ax,syc(indy).^0.5*q,'Normalization','Probability');
legend('X-unc','Y-Unc');
xlabel('Uncertainty in um')
ylabel('Probability');

tN = uitab(tg4,'Title','Number of Photons');
ax = axes(tN);
histogram(ax,N(indy),'Normalization','Probability');
% legend('X-unc','Y-Unc');
xlabel('Photons Detected')
ylabel('Probability');

tllv = uitab(tg4,'Title','LLV');
ax = axes(tllv);
histogram(ax,llv(indy),'Normalization','Probability');
% legend('X-unc','Y-Unc');
xlabel('Log Likelihood Value')
ylabel('Probability');

tiln = uitab(tg4,'Title','LLV./N');
ax = axes(tiln);
histogram(ax,llv(indy)./N(indy),'Normalization','Probability');
% legend('X-unc','Y-Unc');
xlabel('Log Likelihood Value')
ylabel('Probability');

% select data for sx/sy fitting
sxs = sx(indy);
sys = sy(indy);
z0s = z0(indy);
zus = unique(z0s);
for i = 1:numel(zus)
    ind = z0s == zus(i); % select the fits corresponding to current position
    subsx = sxs(ind);
    subsy = sys(ind);
    msx = mean(subsx);
    stx = std(subsx);
    msy = mean(subsy);
    sty = std(subsy);
    ssx(i) = mean(subsx(subsx > msx - 2* stx & subsx < msx + 2* stx));
    ssy(i) = mean(subsy(subsy > msy - 2* sty & subsy < msy + 2* sty));
    
end

% Data Representation
t3 = uitab(tg,'Title','Sigmas');
ax = axes(t3);
plot(ax,z0(indy),sx(indy),'.')
hold on
plot(ax,z0(indy),sy(indy),'.')
% hold off

plot(ax,zus,ssx)
plot(ax,zus,ssy)
hold off

ds = ssx - ssy;

ind = find(abs(ds) == min(abs(ds)));
z0s = zus - zus(ind);
ind = abs(z0s) < 0.7;
z_cal = get_z_params(z0s,ssx(ind),ssy(ind));