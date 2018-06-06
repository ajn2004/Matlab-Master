% Fractals double-exactals
%
%
%
% An updated version of fractal exactals allowing the user to utilize the
% parallelization of the GPU to accelerate rendering and magnification of
% various fractal patterns
% close all;
clearvars; clc

pix = 513;
start = [-1.248,-0.04138];
k = 1;

mag = 1;
wind = 1/mag;
xl = start(1) - wind;
xh = start(1) + wind;
yl = start(2) - wind;
yh = start(2) + wind;

steps = (xh - xl)/pix;

[xgrid, ygrid] = meshgrid(xl:steps:xh,yl:steps:yh);

z0 = complex(xgrid,ygrid);
c = 1;
z = z0;
thresh = 4;
count = z *0;
% for k = 1:500
for  i = 1:500
    z = z.^2 + z0;
    inside = abs(z) <= thresh;
    count = count + inside;
end
imagesc(xgrid(1,:),ygrid(:,1),log(count))
drawnow
M(k) = getframe(gcf);
k = k +1;
