function [im1, um] = make_colorframe(x,y,zmap, gr)
% Make color frames will create a create a frame
[m,n,o] = size(i1); % get image size, this will determine
im1 = i1;
% Because i1 is fixed in size we have to scale the data to match i1
xmax = max(x);
xmin = min(x);
ymax = max(y);
ymin = min(y);

% we know the data set we care about is centered around the origin
% here we calculate the difference in position between max x and max y
dx = xmax - xmin;
dy = xmax - xmin;

if dx <= dy
    um = n/dx; % this is the scale factor to put the data between across 2000 pixels
else
    um = n/dy;
end

% scale the data and remove out of bounds data
x = x*um;
y = y*um;
flag = ( abs(x) > 999.5 | abs(y) > 999.5);

% begin to populate im1
X = -999.5:999.5;
Y = -999.5:999.5;

for i = 1:numel(x)
    if ~flag(i) % only go forward if the point is within bounds
        x_ind = find(X > x(i), 1, 'first');
        y_ind = find(Y > y(i), 1, 'first');
        
       
        % add color value to the pixel identified using RGB convention
        i1(y_ind,x_ind,1) = i1(y_ind,x_ind,1) + zmap(i,1);
        i1(y_ind,x_ind,2) = i1(y_ind,x_ind,2) + zmap(i,2);
        i1(y_ind,x_ind,3) = i1(y_ind,x_ind,3) + zmap(i,3);
    else
    end
end

%% Blurring Kernel setup
ddx = (-50:50)*um;
[x,y] = meshgrid(ddx, ddx);
G = exp(-(x.^2 + y.^2)./(0.5*(gr*um)^2)); % gaussian
gn = G./(sum(G(:))); % normalized

%% Image blurring
for i = 1:3
    im1(:,:,i) = conv2(i1(:,:,i),gn,'same'); % gaussian convolution
end

%% Image Normalization
im1(:,:,1) = im1(:,:,1)./max(max(im1(:,:,1)));
im1(:,:,2) = im1(:,:,2)./max(max(im1(:,:,2)));
im1(:,:,3) = im1(:,:,3)./max(max(im1(:,:,3)));

