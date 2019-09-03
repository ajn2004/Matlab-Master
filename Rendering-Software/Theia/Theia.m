%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theia
%
% Titaness of sight
%
% This is a software suite to build movies of 3D renderings. It will be
% relying heavily on Transformation matrices and the make_xpics software
% group
% AJN 1/27/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars; clc;
mkdir('Render')
% Image settings
dz = 50;
smooth = 1;
grid = 0;
finw = 5; % final width in microns
thresh = 500; % threshold for NN representation
nn_on = 'n';
grayscale = 'y';
%% Movie Resolution Settings
% % 720p movie settings
pixw = 1280;
pixh = 720;

% 1080p movie settings
% pixw = 1920;
% pixh = 1080;

% 4K movie settings
% pixw = 4096;
% pixh = 2160;

%% Zoom Settings
zm_s = 2; % seconds of zoom
fps = 21; % frames per second rate

%% Rotation Settings
tilt = pi/6;  % Tilt around x axis before rotation around y
rot1 = 10;     % seconds for first rotation
rot2 = 4;     % seconds for secont rotation

%% UM bar settings
bx = 50;
by = 50;
um_height = 1/(5*1.61803398875);
vid = VideoWriter('theia_film.avi','Uncompressed AVI');
[fname, fpath] = uigetfile('*tol*');
cd(fpath);

load(fname, 'xf_fixed','yf_fixed','ncoords','q','crlbs','framenumber');
zf = func_build_ramp(ncoords(:,3)*q,framenumber,-20,20,197);
try
xf_all = xf_fixed;
yf_all = yf_fixed;
zf_all = zf/q;
catch lsterr
    load(fname, 'xf_all','yf_all','zf_all','q','crlbs');
end
fin_im = zeros(pixh,pixw,3);
% the basic thought is that we are going to zoom in and rotate around various
% regions in the image and record the path.The image size will stay
% constant so we will be moving more in schrodinger pictures

% preallocate the image space;
% im1 = zeros(pix,pix,3); % create a square canvas for RGB color
finh = pixh*finw/pixw;
% generate color of points
%% Determining Color Information
maxz = ceil(max(zf_all)*q*1000);
minz = floor(min(zf_all)*q*1000);
zs = ceil((maxz-minz)/dz);
z_ind = minz:dz:maxz;
zmap = colormap(jet(zs+1));  % zmap is now a [mol x3] matrix containing color info for each point
% q = q*1000;
%% color selection
for i = 1:numel(zf_all)
    mz = (z_ind - 1000*q*zf_all(i)).^2;  % squared difference between color map and molecule
    zind = find(mz == min(mz),1); % minimum value is effective the z bin
    zc(i,:) = zmap(zind,:);
end

%% Color for NN if it's on
if strcmp(nn_on,'y') || strcmp(nn_on, 'y')
    nn = xf_all*0;
    indy = [];
    for i = 1:numel(xf_all)
        disty = q*1000*((xf_all - xf_all(i)).^2 + (yf_all-yf_all(i)).^2).^0.5;
        nn(i) = min(disty(disty>0));
        if nn(i) < thresh
            indy = [indy;i];
        end
    end
    
    %% Determining Color Information for each molecule
    maxz = ceil(max(nn(indy)));
    minz = floor(min(nn(indy)));
    zs = ceil((maxz-minz)/dz);
    z_ind = minz:dz:maxz;
    zmap = colormap(jet(zs+1));  % zmap is now a [mol x3] matrix containing color info for each point
    
    for i = 1:numel(indy)
        % color selection
        mz = (z_ind - nn(i)).^2;  % squared difference between color map and molecule
        zind = find(mz == min(mz),1); % minimum value is effective the z bin
        zc(i,:) = zmap(zind,:);
    end
    [im1, umb, rad] = make_nnpics(fname, thresh, dz, smooth, grid);
else
%     q = q/1000;
    [im1, umb, rad] = make_3dpics(fname, dz, smooth, grid); % 3D image function
end
% imagesc(im1)

plot(xf_all*q*umb,yf_all*q*umb,'.');
axis image
[x1,y1] = ginput(1); %user selects 2 points as a final frame of the image
[m,n,o] = size(im1);
% pad im1
% im1 = [zeros(m,500,o),im1,zeros(m,500,o)];
% [m,n,o] = size(im1);
% im1 = [zeros(500,n,o);im1;zeros(500,n,o)];
%% BEGIN THE ZOOM

x(1) = (x1) - umb*finw/2;
x(2) = (x1) + umb*finw/2;
y(1) = y1 - umb*finh/2;
y(2) = y1 + umb*finh/2;
nfms = round(zm_s*fps);
lspan = 1:floor(x(1)/nfms):round(x(1));
rspan = x(2):floor((n-x(2))/nfms):n;
dspan = 1:floor(y(1)/nfms):round(y(1));
uspan = y(2):floor((m-y(2))/nfms):m;
vid.FrameRate = fps;
% vid.CompressionRatio = 10;
open(vid);
for k = 1:nfms
    k
    isub = [];
    isub = im1(dspan(k):uspan(nfms-k+1),lspan(k):rspan(nfms-k+1),:);
    [m,n,o] = size(isub);
    fin_im = fin_im*0;
    % Create the scale bar
    um = round(umb);
    
    if um > n/5
        um = um * 0.5;
    end
    
    
    % Resize image for standard sizing
    [m,n,o] = size(isub);
    dh = pixh/m;
    dw = pixw/n;
    
    if m*dw > pixh % check if vert if we scale horz
        dc = dh;
        pixs = pixh;
    else
        dc = dw;
        pixs = pixw;
    end
    
    
    
    
    B = imresize(isub,dc); % resize image into B
    [m,n,o] = size(B);
    % Image Normalization
    B = B - min(B(:));
    B(:,:,1) = B(:,:,1)./max(max(B(:,:,1)));
    B(:,:,2) = B(:,:,2)./max(max(B(:,:,2)));
    B(:,:,3) = B(:,:,3)./max(max(B(:,:,3)));
    bumpm = 1;
    bumpn = 1;
    
    if m > pixh
        m = pixh;
    elseif m < pixh
        bumpm = floor((pixh - m -1)/2)+1;
    end
    
    if n > pixw
        n = pixw;
    elseif n < pixw
        bumpn = floor((pixw - n -1)/2)+1;
    end
    
    fin_im(bumpm:bumpm+m-1,bumpn:bumpn + n -1,:) = B(1:m,1:n,:);
    
    
    um1 = round(um*dc);% Draw the scale bar
    mu1 = round(um1*um_height);
    
    
    fin_im(by:by+mu1-1,bx:bx+um1-1,:) = ones(mu1,um1,3);
    %     for i = by:by+mu1 % spanning y space
    %         for j = bx:bx+um1 % spanning x space
    %             fin_im(i,j,1) = 1;
    %             fin_im(i,j,2) = 1;  % set each channel to 1 causing white
    %             fin_im(i,j,3) = 1;
    %         end
    %     end
    if strcmp(grayscale,'y') || strcmp(grayscale,'Y')
    mfin = max(fin_im,[],3);
    for i = 1:3
        fin_im(:,:,i) = mfin;
    end
    

    end
        imagesc(fin_im(:,:,:));
    
    writeVideo(vid,fin_im);    
    % axis image
    drawnow;
    %     close all
end
% close(vid)
% save('theia_first.mat','fin_im');

disp('pausing');
%% Begin the rotation
% from here on out we'll be dealing with the mathematical coordinates. I am
% making the choice to work entirely in micron space.
% fin_im = [];
% fin_im = zeros(pixs,pixs,3);
% Build micron coordinates
% q = q/1000;
if strcmp(nn_on,'y') || strcmp(nn_on, 'y')
xf_um = xf_all(indy)*q;
yf_um = yf_all(indy)*q;
zf_um = zf_all(indy)*q;
zcc = zc(indy,:);
else
xf_um = xf_all*q;
yf_um = yf_all*q;
zf_um = zf_all*q;
zcc = zc;
end

% select out only coords inside the frame
ind = find( xf_um > x(1)/(umb) & xf_um < x(2)/(umb) & yf_um > y(1)/(umb) & yf_um < y(2)/(umb));
xf = xf_um(ind);
yf = yf_um(ind);
zf = zf_um(ind);
zs = zcc(ind,:);
coords = [xf - mean(xf),yf - mean(yf),zf - mean(zf)]; % centers the coordinates around the origin
% This allows use of
% rotation matrix
% operations

% At this point we can tilt along x and then spin along y to get a
% visualization of what this looks like
nfms = rot1*fps;
theta = 0:2*pi/nfms:2*pi+tilt;


for u = 1:numel(theta)
    tic
    [Rx, Ry, Rz] = rot_mat(theta(u));
    nc = Rx*coords.';
    fin_im = fin_im *0;
    %     plot(nc(1,:),nc(2,:),'.');
    [im2] = make_thisframe(nc(1,:).',nc(2,:).', zs, rad, smooth);
    
    
    
    % Resize image for standard sizing
    [m,n,o] = size(im2);
    dh = pixh/m;
    dw = pixw/n;
    
    if m*dw > pixh % check if vert if we scale horz
        dc = dh;
        pixs = pixh;
    else
        dc = dw;
        pixs = pixw;
    end
    
    
    
    
    B = imresize(im2,dc); % resize image into B
    [m,n,o] = size(B);
    % Image Normalization
    B = B - min(B(:));
    B(:,:,1) = B(:,:,1)./max(max(B(:,:,1)));
    B(:,:,2) = B(:,:,2)./max(max(B(:,:,2)));
    B(:,:,3) = B(:,:,3)./max(max(B(:,:,3)));
    bumpm = 1;
    bumpn = 1;
    
    if m > pixh
        m = pixh;
    elseif m < pixh
        bumpm = floor((pixh - m -1)/2)+1;
    end
    
    if n > pixw
        n = pixw;
    elseif n < pixw
        bumpn = floor((pixw - n -1)/2)+1;
    end
    
    fin_im(bumpm:bumpm+m-1,bumpn:bumpn + n -1,:) = B(1:m,1:n,:);
    
    
    um1 = round(um*dc);% Draw the scale bar
    mu1 = round(um1*um_height);
    
    
    fin_im(by:by+mu1-1,bx:bx+um1-1,:) = ones(mu1,um1,3);
    %     for i = by:by+mu1 % spanning y space
    %         for j = bx:bx+um1 % spanning x space
    %             fin_im(i,j,1) = 1;
    %             fin_im(i,j,2) = 1;  % set each channel to 1 causing white
    %             fin_im(i,j,3) = 1;
    %         end
    %     end
    
    % imagesc(fin_im(:,:,:));
    if strcmp(grayscale,'y') || strcmp(grayscale,'Y')
    mfin = max(fin_im,[],3);
    for i = 1:3
        fin_im(:,:,i) = mfin;
    end
    

    end
        imagesc(fin_im(:,:,:));
    
    writeVideo(vid,fin_im);    
    % axis image
    drawnow;
    % axis image
    % drawnow;
    %     close all
    t(u) = toc;
    ajn_wait(t,u,2*numel(theta));
end
% save('theia_second.mat','fin_im');

clear tic toc t

for u = 1:numel(theta)
    tic
    [Rx, Ry, Rz] = rot_mat(theta(u));
    [Tx, Ty, Tz] = rot_mat(tilt);
    nc = Ry*Tx*coords.';
    fin_im = fin_im*0;
    %     plot(nc(1,:),nc(2,:),'.');
    [im2] = make_thisframe(nc(1,:).',nc(2,:).', zs, rad, smooth);
    
    
    
    % Resize image for standard sizing
    [m,n,o] = size(im2);
    dh = pixh/m;
    dw = pixw/n;
    
    if m*dw > pixh % check if vert if we scale horz
        dc = dh;
        pixs = pixh;
    else
        dc = dw;
        pixs = pixw;
    end
    
    
    
    
    B = imresize(im2,dc); % resize image into B
    [m,n,o] = size(B);
    % Image Normalization
    B = B - min(B(:));
    B(:,:,1) = B(:,:,1)./max(max(B(:,:,1)));
    B(:,:,2) = B(:,:,2)./max(max(B(:,:,2)));
    B(:,:,3) = B(:,:,3)./max(max(B(:,:,3)));
    bumpm = 1;
    bumpn = 1;
    
    if m > pixh
        m = pixh;
    elseif m < pixh
        bumpm = floor((pixh - m -1)/2)+1;
    end
    
    if n > pixw
        n = pixw;
    elseif n < pixw
        bumpn = floor((pixw - n -1)/2)+1;
    end
    
    fin_im(bumpm:bumpm+m-1,bumpn:bumpn + n -1,:) = B(1:m,1:n,:);
    
    um1 = round(um*dc);% Draw the scale bar
    mu1 = round(um1*um_height);
    
    
    fin_im(by:by+mu1-1,bx:bx+um1-1,:) = ones(mu1,um1,3);
    %     for i = by:by+mu1 % spanning y space
    %         for j = bx:bx+um1 % spanning x space
    %             fin_im(i,j,1) = 1;
    %             fin_im(i,j,2) = 1;  % set each channel to 1 causing white
    %             fin_im(i,j,3) = 1;
    %         end
    %     end
%     imagesc(fin_im(:,:,:));
    if strcmp(grayscale,'y') || strcmp(grayscale,'Y')
    mfin = max(fin_im,[],3);
    for i = 1:3
        fin_im(:,:,i) = mfin;
    end
    

    end
%         imagesc(fin_im(:,:,:));
    
    writeVideo(vid,fin_im);    
    % axis image
    drawnow;
    % axis image
    % drawnow;
    %     close all
    t(u) = toc;
    ajn_wait(t,u,numel(theta));
end
close(vid)
% fin_im =[];

% [Rx, Ry, Rz] = rot_mat(pi/2);
% nc = Ry*coords.';
%
% %     plot(nc(1,:),nc(2,:),'.');
%     [im2] = make_thisframe(nc(1,:).',nc(2,:).', zs, rad, smooth);
%     imagesc(im2*40)
%     axis image
%     drawnow

