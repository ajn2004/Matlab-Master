%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center Selector
%   This program renders a density plot of localization data then allows
%   the user to choose centers which will be used for a radial distribution
%   calulation later
%   AJN 2/29/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

pix_size = 20; % Final pixel size in nanometers
p = pix_size/1000;
gauss_std = 1.25; % width of gaussian to convolve with in pixels


[fname, fpath] = uigetfile('*tol.mat');
cd(fpath)

figure('units','normalized','outerposition',[0 0 1 1]);
finfo = dir('*tol.mat');
ind = randperm(numel(finfo));
for lk = 1:numel(finfo)
load(finfo(ind(lk)).name);

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
w = 1;
imagesc(im1);
pl = gca;
hold on
while true
    clearvars points;
    % selects point clicked in plot
    title('Select a new center, press enter to quit')
    bc = waitforbuttonpress;
    if bc == 1
        break      
    else
        points = get(pl,'currentpoint');  
        % assigns selected points to array
        cents(w,1) = points(1,1); % save value of center in um
        cents(w,2) = points(1,2); % save value of center in um
        plot(cents(:,1),cents(:,2),'.r','MarkerSize',35)
        w = w+1;
    end
end
hold off
save_name = [base_name(1:end-4),'_prerdf.mat'];
grid_size = pix_size;
save([fpath,save_name],'im1','cents','xf_all','yf_all','q','p')
clearvars im1 cents xf_all yf_all q
% close all
end