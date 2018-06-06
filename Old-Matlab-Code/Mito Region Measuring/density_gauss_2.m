%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Density Render_2
% AJN 1/26/16
%
% Creates a Density rendering, then convolves the image with a gaussian
% This version contains the gusto needed to finish the mitochondria job
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

pix_size = 20; % Final pixel size in nanometers
p = pix_size/1000;
gauss_std = 1.25; % width of gaussian to convolve with in pixels
cut_off = 8;  % initial density cutoff
ranged = 50; % Range of radii to check for blurring and edge filtering
ranged2 = 5; % range of distances for clustering
overid = 'y';   %use this value to overide chocies
% usev = [ 50 , 5.5];
area_thresh = 1; % Threshold on area selection after identification
type = 'Canny'; % Type of edge filter to be used

% Get and load file
[fname, fpath] = uigetfile('*tol.mat');
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
figure('units','normalized','outerposition', [0 0 1 1])
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

% se = strel('disk',4);
% im3 = imdilate(im2,se);
% im3 = imerode(im3,se);
edge_im = func_gauss_blur(im2,ranged,type);


if ~strcmp(overid,'y');
while true
    radius = input('What radius would you like to use for further processing?');
    
    edge_im = func_gauss_blur(im2,radius,type);
    imagesc(edge_im*(20+mean(im2(:))) + im2);
    
    s = input('Is this a good radius? ','s');
    if strcmp(s,'Y') || strcmp(s,'y')
        break
    end
end
end
p = gca;
while true
    imagesc(edge_im*(20+mean(im2(:))) + im2);
    clearvars points b ys m xs;
    % selects point clicked in plot
    title('Select a new edge, press enter to quit')
    bc = waitforbuttonpress;
    if bc == 1
        break
        
    else
        try
        points = get(p,'currentpoint');
        
        
        % assigns selected points to array
        xs(1) = points(1,1);
        ys(1) = points(1,2);
        
        title('select a second point')
        waitforbuttonpress;
        points = get(p,'currentpoint');
        xs(2) = points(1,1);
        ys(2) = points(1,2);
        
        if xs(2) < xs(1)
            tx = xs(1);
            xs(1) = xs(2);
            xs(2) = tx;
            ty = ys(1);
            ys(1) = ys(2);
            ys(2) = ty;
            clear ty tx
        end
        m = (ys(2)-ys(1))/(xs(2)-xs(1));
        b = ys(2)-m*xs(2);
        xl = xs(1):0.01:xs(2);
        yl = m*xl+b;
        
        index = sub2ind(size(edge_im),round(yl),round(xl));
        edge_im(index) = 1;
        % draws polygon on plot
        catch lsterr
        end
        
    end
end

se = strel('disk',3);
edge_cont = imerode(imdilate(edge_im,se),se);

T = clust_loop(edge_cont,ranged2,overid);
[m, n] = size(edge_im);
[X,Y] = meshgrid(1:n,1:m);
mitos = edge_cont.*0;
% se = strel('disk',4);
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
% imagesc(mitos)
mitos = imopen(mitos,se);
BW = mitos.*im1;
s =regionprops(mitos, 'all');
Centroids = cat(1,s.Centroid);      %
Areas = cat(1, s.Area);             % these lines organize structure data
Diams = cat(1, s.EquivDiameter);    %

%display results of threshold and begin user selection process
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(BW.*im1);
axis equal
title('Circled regions are being selected for with current threshold');
area_ind = find(Areas > area_thresh);
some_cents = Centroids(area_ind,:);
hold on
plot(some_cents(:,1), some_cents(:,2), 'r.');
radii = Diams(area_ind)./2;
circles(some_cents(:,1),some_cents(:,2), radii , 'r');
% hold off
w = numel(some_cents(:,1))+1;
p = gca;
more_cents = some_cents;
more_rads = radii;
while true
    clearvars points;
    % selects point clicked in plot
    title('Select a new center, press enter to quit')
    bc = waitforbuttonpress;
    if bc == 1
        break
        
    else
        
        points = get(p,'currentpoint');
        
        
        % assigns selected points to array
        more_cents(w,1) = points(1,1);
        more_cents(w,2) = points(1,2);
        
        title('select and edge')
        waitforbuttonpress;
        points = get(p,'currentpoint');
        more_rads(w) = ((more_cents(w,1) - points(1,1))^2 + (more_cents(w,2)-points(1,2))^2)^0.5;
        
        % draws polygon on plot
        circles(more_cents(w,1),more_cents(w,2),more_rads(w),'r')
        w = w+1;
    end
end
hold off

save_name = [base_name(1:end-4),'_regions.mat'];
grid_size = pix_size;
save([fpath,save_name],'s','im1','mitos','grid_size','area_ind','more_cents','more_rads','xf_all','yf_all','q')