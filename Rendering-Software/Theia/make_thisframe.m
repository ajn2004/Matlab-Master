function [im1, um] = make_thisframe(xf_all, yf_all, zc, rad, sm)
% Make this frame will create a frame from provided x, y, data

grid_size = 0.5*rad/1000; % pixwidth in um/pix
um = 1/grid_size; % number of pixels required to make a scale bar of 1 um

%% Intializing image grid
xmax = grid_size*ceil(max(xf_all)/grid_size); % Determine maximum x pixels
xmin = grid_size*floor(min(xf_all)/grid_size);
ymax = grid_size*ceil(max(yf_all)/grid_size); % Determine maximum y pixels
ymin = grid_size*floor(min(yf_all)/grid_size);
[Xgrid, Ygrid] = meshgrid(xmin:grid_size:xmax, ymin:grid_size:ymax); % create x and y meshgrid
[m,n] = size(Xgrid);
i1 = zeros(m,n,3); %intialize image variable

%% Determining Color Information
% Color information is predetermined by zc and indexed w/ x and y

%% Populating image grid
% Color will be selected from a pallet and added to an rgb image as
% described above

for i = 1:numel(xf_all) % loop over each molecule
    x_ind = find(Xgrid(1,:) > xf_all(i), 1, 'first') - 1; % determine x bin
    y_ind = find(Ygrid(:,1) > yf_all(i), 1, 'first') - 1; % determine y bin
    
   
    
    % add color value to the pixel identified using RGB convention
    i1(y_ind,x_ind,1) = i1(y_ind,x_ind,1) + zc(i,1);
    i1(y_ind,x_ind,2) = i1(y_ind,x_ind,2) + zc(i,2);
    i1(y_ind,x_ind,3) = i1(y_ind,x_ind,3) + zc(i,3);

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