% Basdi drift preparation
% We want 3-D drift correction, given the size of our data sets we're going
% to choose a subset, or a couple of subsets to correct the drift over, the
% final drift will be an average of the measured drifts in these subsets.
% The drift / frame will be determined by a spline fitting of the final
% drifts.
clearvars; close all; clc;
% % 
% % Load data
[fname, fpath] = uigetfile('*mat');
cd(fpath);
load([fpath, fname],'xf_all','yf_all','zf_all','q','framenum_all');
close all

frames = 5000; % time divisions to find drift over
d_pix = 0.01; % drift pixel size in um
bins = 4;   % number of bins to measure drift over
lat_bin = 2;  % select size of pixels to bin data in um in a lateral projection
               % we'll use this value to create a histogram of the data and
               % select a user defined number of densest locations to
               % perform the drift correction over
% Data conversion to micron space
yf = yf_all*q;
xf = xf_all*q;
zf = (zf_all-min(zf_all))*q;
               
% Select subset of data
mx = ceil(max(xf)/lat_bin); % get max values
my = ceil(max(yf)/lat_bin);
mf = max(framenum_all);
% convert into lat_bin_pixels
yb = yf/lat_bin;
xb = xf/lat_bin;
 % Create a 2D histogram to determine density bins
% [xg, yg] = meshgrid(-0:lat_bin:mx,-0:lat_bin:my);
I = zeros(my+2,mx+2);
for i = 1:numel(xb)
%     [col] = find( abs(xg(1,:) - xb(i)) < 0.5, 1,'First') ;
%     [row] = find( abs(yg(:,1) - yb(i)) < 0.5, 1,'First') ;

    I(round(yb(i))+1,round(xb(i))+1) = I(round(yb(i)+1),round(xb(i))+1) + 1;

end
  
% select the densest bins
IS = sort(I(:)); % sort the image from smallest to largest
tbs = IS(end-bins+1:end);  % densest bins value

% Collect index of brightest bins
row = [];
col = [];
for i =1 :numel(tbs)
    [row(i),col(i)] = find(I == tbs(i),1);
end

% grab the coords of the densest pixels based on index
C = {};
imagesc(I)
hold on
for i = 1:numel(row)
    ind = abs(xb - col(i)+1) < 0.5 & abs(yb - row(i)+1)< 0.5; % create index of 
    C{i} = [xf(ind),yf(ind),zf(ind),framenum_all(ind)];
    plot(xb(ind)+1,yb(ind)+1,'.');
%         plot3(xf(ind),yf(ind),zf(ind),'.');
%     hold on
%     axis equal
end
clear row col I ind tbs IS
% at this point we have a cell array of coords relating to the 'bins' densest
% areas of the image, we will perform the drift correction over these bins
% C now contains all variables of interest for correction

% Drift Estimation
% In this section we're going to attempt to estimate the Drift of the
% regions chosen
% for i = 1:numel(C)
frm_dvs = ceil(mf/frames);
for i = 1:bins % start with the first case then we'll loop over all C
    data = [];
    data = C{i}; % C data is stored as [x, y, z, frame]
    xs = data(:,1);
    ys = data(:,2);
    zs = data(:,3);
    fs = data(:,4);
    % heigh and width section
    wxy = round((ceil(max(xs))+1)/d_pix); % for x-y
    hxy = round((ceil(max(ys))+1)/d_pix);
    wxz = round((ceil(max(xs))+1)/d_pix); % for x-z
    hxz = round((ceil(max(zs))+1)/d_pix);
    wyz = round((ceil(max(ys))+1)/d_pix); % for y-z
    hyz = round((ceil(max(zs))+1)/d_pix);
    % create array of cells for basdi drift
    Cxy = {};
    Cxz = {};
    Cyz = {};    
    
    % build cell data
    for j = 1:frm_dvs
        ind = data(:,4) >= (j-1)*frames & data(:,4) <= j*frames;
        if ~isempty(ind)
            Cxy{j} = [round(ys(ind)/d_pix),round(xs(ind)/d_pix)];
            Cxz{j} = [round(zs(ind)/d_pix),round(xs(ind)/d_pix)];
            Cyz{j} = [round(zs(ind)/d_pix),round(ys(ind)/d_pix)];
        else
            Cxy{j} = [];
            Cxz{j} = [];
            Cyz{j} = [];
        end
    end
% 	sxy(i) = BaSDI_main(Cxy,hxy,wxy);
%     sxz(i) = BaSDI_main(Cxz,hxz,wxz);
%     syz(i) = BaSDI_main(Cyz,hyz,wyz);
        
end