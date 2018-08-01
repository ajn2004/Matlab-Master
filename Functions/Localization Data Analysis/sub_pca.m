function [rx, ry, ind] = sub_pca(xf,yf)
plot(xf,yf,'.')
[x,y] = ginput(2);
ind = xf < max(x) & xf > min(x);
ind = ind & yf < max(y) & yf > min(y);
xs = xf(ind) - mean(xf(ind));
ys = yf(ind) - mean(yf(ind));
[M,V] = eig(cov(xs,ys));
rots = M*[xs.';ys.'];
rx = rots(1,:);
ry = rots(2,:);
