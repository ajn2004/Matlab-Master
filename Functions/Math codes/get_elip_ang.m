function [top_ang, C] = get_elip_ang(im1,sigx,sigy)
% find the best angle to fit an eliptical gaussian to im1
% im1 must be a square for this to work
offs = min(im1(:)); % approx offset is the minimum pixel value
peak = max(im1(:)); % approx peak value is maximum pixel value
top_ang = 0; % preallocate top angle
mC = 10000000000000000; % set min cost variable to absurd high
[m] = size(im1); % get size of image

% handle cases if m is even or odd
if m/2 == round(m/2) % if m is even m/2 is an integer and = round(m/2)
    n = m/2;
    [X,Y] = meshgrid(-n+0.5:n-0.5);
else % if the above doesn't hold m must be odd
    n = (m-1)/2;
    [X,Y] = meshgrid(-n:n);
end
xcm = sum(sum(im1.*X./sum(im1(:))));
ycm = sum(sum(im1.*Y./sum(im1(:))));
for i = 0:179 % loop around in degrees
    rx = (X-xcm)*cos(deg2rad(i)) - (Y-ycm)*sin(deg2rad(i)); % rotate underlying grid
    ry = (X-xcm)*sin(deg2rad(i)) + (Y-ycm)*cos(deg2rad(i));
    imx = peak*exp(-rx.^2/(2*sigx^2) - ry.^2/(2*sigy^2)) + offs; % approximate gaussian
    C(i+1) = sum(sum((imx - im1).^2)); % determin cost as sum of squares
    if C(i+1) < mC % check for minimum cost
        top_ang = deg2rad(i);
        mC = C(i+1);
    end
end
    
    