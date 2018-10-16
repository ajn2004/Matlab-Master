% Attempting to figure out the tilt profile
xt = [];
yt = [];
zt = [];
for i = 1:max(pfs)
    ind = pfs == i;
    xt = [xt;xfs(ind)-mean(xfs(ind))];
    yt = [yt;yfs(ind)-mean(yfs(ind))];
    zt = [zt;zfs(ind)];
end