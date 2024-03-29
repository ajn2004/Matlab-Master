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
    A{i} = rollingball(readtiff(files(i).name),5,4);                                         % store total images into variable
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
psf = [];
% gpu fit of PSFs
% NOTE TO USER: because we assume an axis of elipticity that
% need not be parallel to the axis of the camera we rotate our
% grid to perform the fit. This rotation is corrected for on
% the back end and the X-Y positions returned are understood to
% be on the axis of the camera. However, the sigma values exist
% along the rotated or prime axis. This is representated in the
% following comments once and is to be understood that
% Sigma-x/y corresponds to the prime axis and not the camera
% axis
fms = 1:o; % input framenum
for i = 1:numel(psfs)                                                       % Loop over all identified ROIs
    [fits, crlb, lv,fnout] = slim_locs(psfs{i},fms,zeros(o,2),ang,50,100);  % Perform Fit
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
    psf = [psf; fits(:,1)*i./fits(:,1)];                                            % Log the PSF the fit is associated with
    llv = [llv;-abs(lv)];                                                   % Log Likelihood Value
%     fnum = [fnum;fnout - disp(i)];                                          % Correct the Frame number based off correlation result
    fnum = [fnum;fnout]; 
end
z0 = fnum*step/1000;                                                        % Populate Z-positions
% Tolerance Step
%(User can adjust these values to change fidelity of data used to determine calibration results)
% For use with 10-3-18 scans
% indy = llv./N > -0.1 & N >0 & N < 10000 & abs(syc) < 0.01 & abs(sxc) < 0.01 ;
indy = llv./N > -0.15 & N >0 & N < 100000 & abs(syc).^0.5*q < 0.006 & abs(sxc).^0.5*q < 0.006 ;

%% Data Representation

% Display X-Y Position
t4 = uitab(tg,'Title','Fitting Outputs');
tg4 = uitabgroup(t4);
txy = uitab(tg4,'Title','X-Y Position');
ax = axes(txy);
plot(ax,xf(indy),yf(indy),'.');
axis equal
xlabel('X-Fit')
ylabel('Y-Fit');

% Display X-Y Uncertainties
tlp = uitab(tg4,'Title','X-Y Uncertainty');
ax = axes(tlp);
histogram(ax,xfc(indy).^0.5*q,'Normalization','Probability');
hold on
histogram(ax,yfc(indy).^0.5*q,'Normalization','Probability');
legend('X-unc','Y-Unc');
xlabel('Uncertainty in um')
ylabel('Probability');
hold off

% Display Uncertainties in Sigma Values
tsxy = uitab(tg4,'Title','Sigma Uncertainties');
ax = axes(tsxy);
histogram(ax,sxc(indy).^0.5*q,'Normalization','Probability');
hold on
histogram(ax,syc(indy).^0.5*q,'Normalization','Probability');
legend('X-unc','Y-Unc');
xlabel('Uncertainty in um')
ylabel('Probability');

% Display Histogram of Number of Photons
tN = uitab(tg4,'Title','Number of Photons');
ax = axes(tN);
histogram(ax,N(indy),'Normalization','Probability');
xlabel('Photons Detected')
ylabel('Probability');

% Display Histogram of Log-Likelihood Values
tllv = uitab(tg4,'Title','LLV');
ax = axes(tllv);
histogram(ax,llv(indy),'Normalization','Probability');
xlabel('Log Likelihood Value')
ylabel('Probability');

% Display 'goodness of fit' metric LLV/N
tiln = uitab(tg4,'Title','LLV./N');
ax = axes(tiln);
histogram(ax,llv(indy)./N(indy),'Normalization','Probability');
xlabel('Log Likelihood Value')
ylabel('Probability');

%% Z-Calibration
% select subset of data that passed tolerances for sx/sy fitting
sxs = sx(indy);                                                             % Sigma X
sys = sy(indy);                                                             % Sigma Y
z0s = z0(indy);                                                             % Z-values
zus = unique(z0s);                                                          % Create a Unique list of Z-positions for indexing purposes
for i = 1:numel(zus)                                                        % Loop over unique Z-Positions
    ind = z0s == zus(i);                                                    % select the fits corresponding to a particular Z-Position
    subsx = sxs(ind);                                                       % Create subset of Sigma-x corresponding to this Z-pos
    subsy = sys(ind);                                                       % Create subset of Sigma-y corresponding to this Z-pos
    msx = mean(subsx);                                                      % Determine Average Value of Sigma-x at this Z-pos
    stx = std(subsx);                                                       % Determine standard deviation of Sigma-x at this Z-pos
    msy = mean(subsy);                                                      % Determine Average Value of Sigma-y at this Z-pos
    sty = std(subsy);                                                       % Determine standard deviation of Sigma-x at this Z-pos
    ssx(i) = mean(subsx(subsx >= msx - 2* stx & subsx <= msx + 2* stx));      % Select only values within 2 standard deviations of the average for X'
    ssy(i) = mean(subsy(subsy >= msy - 2* sty & subsy <= msy + 2* sty));      % Select only values within 2 standard deviations of the average for Y'
    
end

% Find plane of least confusion
ds = ssx - ssy;                                                             % Subtract sigma-y from sigma-x
ind = find(abs(ds(40:end-40)) == min(abs(ds(40:end-40))));                                        % Find minimum of absolute value
z0s = zus - zus(ind(1)+40);                                                       % Call the found index the 0 point

% Data Representation of Sigma and Z space
t3 = uitab(tg,'Title','Sigmas');
ax = axes(t3);
plot(ax,z0(indy)-zus(ind(1)+40),sx(indy),'.')                                  % Plot Toleranced Sigma-x data
hold on
plot(ax,z0(indy)-zus(ind(1)+40),sy(indy),'.')                                  % Plot Toleranced Sigma-y data
plot(ax,z0s,ssx)                                                            % Plot average sigma-x
plot(ax,z0s,ssy)                                                            % Plot average sigma-y

% Determining Z parameters
ind = abs(z0s) < 0.7; % Limit
z_cal = get_z_params(z0s(ind),ssx(ind),ssy(ind));

% Display results of Z-calibration
yx = z_cal_fit(z0s(ind),z_cal(1:5));    % Determine fitted Sig-x Values
yy = z_cal_fit(z0s(ind),z_cal(6:end));  % Determine fitted Sig-y values
% Overlay result on scatted / average sig-x plot
plot(ax,z0s(ind),yx,'gx')
plot(ax,z0s(ind),yy,'gx')
hold off

%% Determine Necessary x-y-z correction
% It's known that in astigmatism there may be a slight slant that
% exists over the Z direction, in this section we'll address that

zf_um = getdz(sx,sy,z_cal); % Get Z-values
indy = indy & abs(zf_um) <0.6;
d3 = uitab(tg4,'Title','3-D Positions');
ax = axes(d3);
scatter3(xf(indy),yf(indy),zf_um(indy)/q,[],psf(indy));
axis equal
xlabel('Lateral-X')
ylabel('Lateral-Y');
zlabel('Axial-Z');

% Display z0 true v zf
tt4 = uitab(tg,'Title','Final Corrections');
tg5 = uitabgroup(tt4);
tf = uitab(tg5,'Title','Z0 v. Zf');
ax = axes(tf);
plot(ax,z0(indy),zf_um(indy),'.')
a = polyfit(z0(indy),zf_um(indy),1);  % The slope of this distribution gives us the 'correction' for absolute Z
hold on
plot(ax,z0(indy),a(1)*z0(indy)+a(2),'g')
legend('Data','Fit','Location','North');
hold off
xlabel('Stage Position');
ylabel('Found Position');

% grab subsets
dist = 0.5;
indy = indy & abs(zf_um) < dist;
xfs = xf(indy);
yfs = yf(indy);
zfs = zf_um(indy)/q;
pfs = psf(indy);
zs = (min(zfs*q):0.02:max(zfs*q))/q;  % zfs is in pixels, are values are in pixels right now

% Align the PSFs by their average positions
xt = [];
yt = [];
zt = [];
for i = 1:max(pfs)
    ind = pfs == i;
    xt = [xt;xfs(ind)-mean(xfs(ind))];
    yt = [yt;yfs(ind)-mean(yfs(ind))];
    zt = [zt;zfs(ind)];
end
% indy = indy & xfc.^0.5*q < 0.005 & yfc.^0.5*q < 0.005;
% indy = indy & abs(xf) < 0.5 & abs(yf) < 0.5;
% Grab average position in the aligned psfs
xsel = [];
ysel = [];
zsel = [];
for i = 1:numel(zs) - 1
    ind1 = zt >= zs(i) & zt <=zs(i+1);
    ind1 = ind1 & xt > -0.2 & xt < 0.5;
    ind1 = ind1 & yt > -0.25 & yt < 0.25;
    xts = xt(ind1);
    yts = yt(ind1);
    zts = zt(ind1);
    mx = mean(xts);
    stx = std(xts);
    my = mean(yts);
    sty = std(yts);
    ind2 = xts >= mx - 1*stx & xts <= mx +1*stx & yts >= my - 1*sty & yts <= my +1*sty;
    xsel = [xsel;xts(ind2)];
    ysel = [ysel;yts(ind2)];
    zsel = [zsel;zts(ind2)];
    xfm(i) = mean(xts(ind2));
    yfm(i) = mean(yts(ind2));
    zfm(i) = mean(zts(ind2));
end

splz = (-dist:0.001:dist)/q;
xtilt = gausssmooth(xfm,5,10);
ytilt = gausssmooth(yfm,5,10);
splx = spline(zs(1:end-1) + mean(diff(zs))/2,xtilt,splz);
sply = spline(zs(1:end-1) + mean(diff(zs))/2,ytilt,splz);

ts = uitab(tg5,'Title','Z-Tilt w/ Spline');
ax = axes(ts);
plot3(ax,xfm,yfm,zfm,'.','Color',[0,1,0],'MarkerSize',20);
hold on
plot3(ax,xt,yt,zt,'.','Color',[0.1,0.1,0.1],'MarkerSize',2);
plot3(ax,xsel,ysel,zsel,'.','Color',[0,0,1],'MarkerSize',8);
plot3(ax,splx,sply,splz,'Color',[1,0,0],'LineWidth',10);
hold off
xlabel('Lat-X');
ylabel('Lat-Y');
zlabel('Axi-Z');
legend('Avg Pts','Spline Curve');

% Make Drift Correction
zf = zf_um/q;
xc = spline(zs(1:end-1) + mean(diff(zs))/2,xtilt,zf);
yc = spline(zs(1:end-1) + mean(diff(zs))/2,ytilt,zf);

xf_c = xf - xc;
yf_c = yf - yc;
zf_c = zf/(a(1));
tc = uitab(tg5,'Title','Corrected localizations');
ax = axes(tc);
scatter3(ax,xf_c(indy),yf_c(indy),zf_c(indy),[],psf(indy));
colormap('jet')
xlabel('Lat-X');
ylabel('Lat-Y');
zlabel('Axi-Z');
axis equal
hold off
cal.z_cal = z_cal;
cal.a = a;
cal.tilt.x = xtilt;
cal.tilt.y = ytilt;
cal.ang = ang;
cal.tilt.zs = zs;
save('z_calib.mat','cal')
