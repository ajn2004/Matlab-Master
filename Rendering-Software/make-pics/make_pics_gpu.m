function i1 = make_pics_gpu(coord,pixs,dim)
% Matlab wrapper for density_hist.cu
% Coords and bin are expected to be in the SAME UNITS!

%% Put all coords in 'positive' space
for i = 1:2
    coord(:,i) = coord(:,i) - min(coord(:,i));
end

xf = coord(:,1);
yf = coord(:,2);

xpix = pixs(2);
ypix = pixs(1);
xl = [min(xf), max(xf)];
yl = [min(yf), max(yf)];
m = pixs(1);
n = pixs(2);
if dim == 2
    scale = n/xl(2);
else
    scale = m/yl(2);
end

xf = ((xf - xl(1))*scale);
yf = ((yf - yl(1))*scale);
xf = xf + 960/2;

m = int64(m);
n = int64(n);


[i1] = make_gpupics(single([xf,yf]),int64([m;n]));




