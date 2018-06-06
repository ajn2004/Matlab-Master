% Fractals double-exactals
%
%
%
% An updated version of fractal exactals allowing the user to utilize the
% parallelization of the GPU to accelerate rendering and magnification of
% various fractal patterns
% close all;
% clearvars; clc

pix = 500;
% % for pix = 1:1000
% [xgrid, ygrid] = meshgrid(-round(pix/2):round(pix/2),-round(pix/2):round(pix/2));
m = 387;
n = 87;
start = [xgrid(m,n),ygrid(m,n)];
% start = [-0.2003, 0.5437];
k = 1;
iter = 10000;
i = 0;
count = 1;
% while true

for i = 0:0.01:18
mag = 10^i;
wind = 1/mag;
xl = start(1) - wind;
xh = start(1) + wind;
yl = start(2) - wind;
yh = start(2) + wind;

steps = (xh - xl)/pix;

[xgrid, ygrid] = meshgrid(xl:steps:xh,yl:steps:yh);

z0 = complex(xgrid,ygrid);
z = z0;
thresh = 4;

im1 = frac_exac(xgrid, ygrid, xgrid, ygrid, thresh, iter);
imagesc(log(im1));
axis image
title(['log(Magnifcation)', num2str(i)]);
drawnow
M(count) = getframe(gcf);
count = count +1;
% [x, y] = ginput(1);
% m = round(y);
% n = round(x);
% start = [xgrid(m,n),ygrid(m,n)];
% i = i+1;
% start = [-0.5, 0];
end
% m = round(numel(xgrid(:,1))/2);
% n = round(numel(xgrid(1,:))/2);
%     end
% end