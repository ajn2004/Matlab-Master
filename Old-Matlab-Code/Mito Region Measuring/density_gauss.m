%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Density Render
% AJN 1/26/16
%
% Creates a Density rendering, then convolves the image with a gaussian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

pix_size = 20; % Final pixel size in nanometers
p = pix_size/1000;
gauss_std = 1.25; % width of gaussian to convolve with in pixels
cut_off = 8;  % initial density cutoff
ranged = 1:50; % Range of radii to check for blurring and edge filtering
ranged2 = 0.5:0.5:10; % range of distances for clustering
type = 'Canny'; % Type of edge filter to be used

% Get and load file
[fname, fpath] = uigetfile('*.mat');
addpath(pwd);
cd(fpath);
load(fname);

% prepare variables for rendering
xf_in = xf_all*q;
yf_in = yf_all*q;

% determine maximum plotted position in microns
max_x = ceil(max(xf_in)/p)*p+p;
max_y = ceil(max(yf_in)/p)*p+p;
[Xgrid, Ygrid] = meshgrid(0:p: max_x,0:p: max_y);
dens = zeros(size(Xgrid));
[m, n] = size(Xgrid);
for i = 1:numel(xf_in)
    x_ind = find(Xgrid(1,:) > xf_in(i), 1, 'first') - 1;
    y_ind = find(Ygrid(:,1) > yf_in(i), 1, 'first') - 1;
    dens(y_ind,x_ind) = dens(y_ind,x_ind) + 1;
end

% imagesc(dens)
[X, Y] = meshgrid(-10:10,-10:10);
gauss = exp(-2*((X).^2 + (Y).^2)./(gauss_std*2)^2);

im1 = conv2(dens,gauss,'same');

%% Density thresholding
while true
    [row, col] = find(im1 < cut_off);
    im2 = im1;
    
    for j = 1:numel(row)  % I don't know why, but without this loop program freezes
        im2(row(j),col(j)) = 0; % loop over indecies and set to 0
    end
    
    
    imagesc(im2);
    title([num2str(cut_off) , 'is the current threshold']);
    drawnow
    s = input('Was this threshold good? ', 's');
    if strcmp(s,'Y') || strcmp(s,'y')
        break
    end
    cut_off = input('Choose a new threshold for density');
end
edge_im = func_gauss_blur(im2,ranged,type);

while true
    radius = input('What radius would you like to use for further processing?');
    
    edge_im = func_gauss_blur(im2,radius,type);
    imagesc(edge_im*(20+mean(im2(:))) + im2);
    
    s = input('Is this a good radius? ','s');
    if strcmp(s,'Y') || strcmp(s,'y')
        break
    end
end


se = strel('disk',3);
edge_cont = imerode(imdilate(edge_im,se),se);

T = clust_loop(edge_cont,ranged2);
[m, n] = size(edge_im);
[X,Y] = meshgrid(1:n,1:m);
mitos = edge_cont.*0;
se = strel('disk',4);
for i = 1:max(T(:,1))
    indpol = find(T(:,1) == i);
    templ = mitos.*0;
    if ~isempty(indpol)
        templ(T(indpol,2),T(indpol,3))=edge_cont(T(indpol,2),T(indpol,3));
        templ = imdilate(templ,se);
        templ = imfill(templ,'holes');
        templ = imerode(templ,se);
%         templ(T(indpol,2),T(indpol,3))=im1(T(indpol,2),T(indpol,3));
        mitos = logical(mitos + templ);
        
    end
end
imagesc(mitos)