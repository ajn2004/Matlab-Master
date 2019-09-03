function [i1, i2, i3] = func_3D_dens(coord,bin,radius)
% Matlab wrapper for density_hist.cu
% Coords and bin are expected to be in the SAME UNITS!
for i = 1:3
    coord(:,i) = coord(:,i) - min(coord(:,i));
end
xf = coord(:,1)*1000;
yf = coord(:,2)*1000;
zf = coord(:,3)*1000;
radius = radius *1000;
sbin = bin;
xl = [min(xf), max(xf)];
yl = [min(yf), max(yf)];
zl = [min(zf), max(zf)];
m = (floor(diff(yl)/sbin)+1);
n = (floor(diff(xl)/sbin)+1);
o = (floor(diff(zl)/sbin)+1);
% i1 = int32(zeros(m,n,o)); % Create 3D histogram
% xf = (round((xf - xl(1))*n/xl(2)));
% yf = (round((yf - yl(1))*m/yl(2)));
% zf = (round((zf - zl(1))*o/zl(2)));
xf = (((xf - xl(1))*n/xl(2)));
yf = (((yf - yl(1))*m/yl(2)));
zf = (((zf - zl(1))*o/zl(2)));
rads = single(radius*n/xl(2));
% Populate the 3D Histogram
% plot3(xf,yf,zf,'.')
% drawnow
% axis equal
m = int64(m);
n = int64(n);
o = int64(o);


% tic
[i1] = density_hist(xf,yf,zf,int64([m;n;o]),rads);
% toc



