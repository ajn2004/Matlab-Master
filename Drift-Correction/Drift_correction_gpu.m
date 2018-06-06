%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Drift Correction on GPU
%
%  This is a script to address drift in FPALM images analyzed with Neural
%  Quhzx and other GPU oriented codes. The idea is that it will follow the
%  same math as the drift correction gui but
%
%
%
%
% AJN 7-21-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

pix_size = 25; % Final pixel size in nanometers
chunk_size = 1000; % number of frames to construct a partial render for correlation
max_frames = 10000;
p = pix_size/1000;
gauss_std = 1.75; % width of gaussian to convolve with in pixels


[fname, fpath] = uigetfile('*tol.mat');
cd(fpath)

% figure('units','normalized','outerposition',[0.5 0 0.5 1]);
% finfo = dir('*tol.mat');
% randin = randperm(numel(finfo));

load(fname)

% prepare variables for rendering
xf_in = xf_all*q;
yf_in = yf_all*q;
% max_frames = max(framenum_all);
chunks = floor((max_frames -1)/chunk_size)+1;
[X, Y] = meshgrid(-10:10,-10:10);
gauss = exp(-2*((X).^2 + (Y).^2)./(gauss_std*2)^2);
% Build the images for subsequent correlations
for j = 1: chunks 
    % build xf and yf part for density plot rendering
    frame1 = (j-1)*chunk_size + 1;
    frame2 = j*chunk_size+1;
    ind = framenum_all >= frame1 & framenum_all <= frame2;
    xf_part = xf_in(ind);
    yf_part = yf_in(ind);
    
    % determine maximum plotted position in microns
    max_x = ceil(max(xf_in)/p)*p+p;
    max_y = ceil(max(yf_in)/p)*p+p;
    [Xgrid, Ygrid] = meshgrid(0:p: max_x,0:p: max_y);
    dens = zeros(size(Xgrid));
    [m, n] = size(Xgrid);
    for i = 1:numel(xf_part)
        x_ind = find(Xgrid(1,:) > xf_part(i), 1, 'first') - 1;
        y_ind = find(Ygrid(:,1) > yf_part(i), 1, 'first') - 1;
        dens(y_ind,x_ind) = dens(y_ind,x_ind) + 1;
    end

    
    dens = conv2(dens,gauss,'same');
    cmd_str = ['im',num2str(j),'=dens;'];
    eval(cmd_str);
%     cmd_str = ['imagesc(im',num2str(j),');'];
%     eval(cmd_str);
%     pause(0.2);
end
[m,n] = size(im1);
% perform cross correlations and find shift
drifts = [0,0];
for i = 1:chunks-1 % start on frame i and correlate with frame i+1
    cmd_str = ['the_coords = xcorr2(im',num2str(i), ', im',num2str(i+1),');'];
    eval(cmd_str);
%     sub_coords = the_coords(m-10:m+10,n-10:n+10);
    [row, col] = find(imgaussfilt(the_coords,2) == max(max(imgaussfilt(the_coords,2))));
    drifts = [drifts; pix_size*(col(1) - n), pix_size*(row(1) - m)];
end

% correct local drift to absolute drift
for i = 2:numel(drifts(:,1))
    drifts(i,:) = drifts(i-1,:) + drifts(i,:);
end

% determining functional form of drift in x and y
%in this section we will loop over a 9th order polynomial 
n = numel(drifts(:,1));
x = (0.5:1:chunks-0.5)*chunk_size;
x = x.';

for i = 0:n-2 % here we are going to go to number of data points - 2 so we are avoid over fitting
    px = polyfit(x,drifts(:,1),i);
    py = polyfit(x,drifts(:,2),i);
    % calculate sum of squared differences to determine best polynomial to
    % use
    srx(i+1) = sum((drifts(:,1) - polyval(px,x)).^2)/(n- i -1);
    sry(i+1) = sum((drifts(:,2) - polyval(py,x)).^2)/(n- i -1);
end

% use var to find what order the polynomial is no longer significantly
% changing
flagx = 0;
for i = 1:numel(srx)
    varx(i) = var(srx(i:end));
    if i>2
        if varx(i)/varx(i-1) > varx(i-1)/varx(i-2)
            orx = i-2;
            flagx = 1;
            break
        end
    end
end
flagy = 0;
for i = 1:numel(srx)
    vary(i) = var(sry(i:end));
    if i>2
        if vary(i)/vary(i-1) > vary(i-1)/vary(i-2)
            ory = i-2;
            flagy = 1;
            break
        end
    end
end
if flagx == 0
    orx = 0;
end
if flagy == 0
    ory = 0;
end

% calculate fit 
px = polyfit(x,drifts(:,1),orx);
py = polyfit(x,drifts(:,2),ory);

% display fit
figure
subplot(2,2,1);
plot(x,drifts(:,1),'.b');
hold on
plot(x,polyval(px,x));
hold off
xlabel('frame chunks'); 
ylabel('drift in nm'); 
title('X drift');

subplot(2,2,3);
plot(x,drifts(:,2),'.b');
hold on
plot(x,polyval(py,x));
xlabel('frame chunks'); 
ylabel('drift in nm'); 
title('Y drift');

subplot(2,2,[2,4]);
plot(drifts(:,1),drifts(:,2));
xlabel('xdrift');
ylabel('ydrift');
title('Total Drift');
axis equal

% correct data
xdrift = polyval(px,framenum_all)/(1000*q);
ydrift = polyval(py,framenum_all)/(1000*q);

xf_fixed = xf_all-xdrift;
yf_fixed = yf_all-ydrift;

plot(xf_all,yf_all,'.r');
hold on
plot(xf_fixed,yf_fixed,'.b');
hold off