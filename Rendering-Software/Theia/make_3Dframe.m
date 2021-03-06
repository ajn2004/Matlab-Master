function [im1, um] = make_3Dframe(fdname, varargin)
% This is a rednering function that will create an rgb([0 1]) image
% Variables includ dz, smoothing, and grid_size
% code adapted from show_pics in da_support folder
load(fdname); 
numvar = length(varargin);
if numvar < 1
    dz = 20;
    sm = 1;
elseif numvar == 1
    dz = varargin{1};
    sm = 1;
elseif numvar >= 2
    dz = varargin{1};
    sm = varargin{2};
end
%% Figuring out pixel size
q = q*1000; % convert q to nm
rad = mean((xf_crlb.*yf_crlb).^0.25)*q; % all distances in this are determined off of loc_unc
grid_size = 0.5*rad; % pixwidth in nm/pix
um = 1000/grid_size; % number of pixels required to make a scale bar of 1 um
if numvar >= 3
    gride = varargin{3};
    if gride > 0
        grid_size = gride;
    end       
end

plot(xf_all,yf_all,'.');
[x,y] = ginput(1);
% Recenter
xf_all = xf_all - x;
yf_all = yf_all - y;
close all
%% Intializing image grid
xmax = grid_size*ceil(max(xf_all)*q/grid_size); % Determine maximum x pixels
ymax = grid_size*ceil(max(yf_all)*q/grid_size); % Determine maximum y pixels
[Xgrid, Ygrid] = meshgrid(0:grid_size:xmax, 0:grid_size:ymax); % create x and y meshgrid
i1 = zeros(round(ymax/grid_size),round(xmax/grid_size),3); %intialize image variable

%% Determining Color Information
maxz = ceil(max(zf_nm)); 
minz = floor(min(zf_nm));

zs = ceil((maxz-minz)/dz);
z_ind = minz:dz:maxz;
zmap = colormap(jet(zs+1));

%% Populating image grid
% Color will be selected from a pallet and added to an rgb image as
% described above
for i = 1:numel(xf_all) % loop over each molecule
    x_ind = find(Xgrid(1,:) > xf_all(i)*q, 1, 'first') - 1; % determine x bin
    y_ind = find(Ygrid(:,1) > yf_all(i)*q, 1, 'first') - 1; % determine y bin
    
    % color selection
    mz = (z_ind - zf_nm(i)).^2;  % squared difference between color map and molecule
    zind = find(mz == min(mz),1); % minimum value is effective the z bin
    
    % add color value to the pixel identified using RGB convention
    i1(y_ind,x_ind,1) = i1(y_ind,x_ind,1) + zmap(zind,1);
    i1(y_ind,x_ind,2) = i1(y_ind,x_ind,2) + zmap(zind,2);
    i1(y_ind,x_ind,3) = i1(y_ind,x_ind,3) + zmap(zind,3);

end

%% Blurring Kernel setup
[x,y] = meshgrid(-50:50,-50:50);
G = exp(-(x.^2 + y.^2)./(0.5*(sm*rad)^2)); % gaussian
gn = G./(sum(G(:))); % normalized

%% Image blurring
for i = 1:3
    im1(:,:,i) = conv2(i1(:,:,i),gn,'same'); % gaussian convolution
end

%% Image Normalization
im1(:,:,1) = im1(:,:,1)./max(max(im1(:,:,1)));
im1(:,:,2) = im1(:,:,2)./max(max(im1(:,:,2)));
im1(:,:,3) = im1(:,:,3)./max(max(im1(:,:,3)));

clear sim1 G gn x y drifts fpath framenum_all i iln iloc imgf llv N N_crlb nfiles off_all
clear off_crlb pixtpho pixw q sigx_all sigx_crlb sigy_all sigy_crlb simgf
clear total_mol total_molecues xdrift xf_all xf_crlb xf_fixed xf_part ydrift yf_all yf_crlb
clear yf_fixed yf_max zf_nm zf_nm_temp
save(['render\',fdname(1:end-4),'_render.mat']);