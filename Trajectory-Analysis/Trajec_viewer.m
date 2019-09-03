% Trajec Viewer
close all
xf_all = ncoords(:,1);
yf_all = ncoords(:,2);
zf_all = ncoords(:,3);
xpo = 0.04;
mmax = [0 10];
% for lk = 1:60
lk = 0;
zf_all = func_shift_correct(ncoords(:,3)*q,framenumber,1)/q;
% zf_all = ncoords(:,3);
zf_all = zf_all(:);
scl = 1;
xf_all = xf_fixed;
yf_all = yf_fixed;
% zf_all = ncoords(:,3);
mN = 2000;
msnr = 10;
% s = scatter3(xf_all*q,yf_all*q,zf_all*q,-scl*llv./fits(:,3));
snr = fits(:,3)./(fits(:,3) + (pixw*2+1)^2*fits(:,6)).^0.5;
% plot3(q*xf_all,q*yf_all,q*zf_all,'.');
axis equal
hold on
dx = [];
dy = [];
dz = [];
n = [];
zd = 0;
ddz = [];
mz = [];
mf = [];

for i = 1:numel(trajec)
    ind = trajec(i).t;
    nums(i) = numel(ind);
    indy = snr(ind) < msnr;
    ind(indy) = [];
    rgb = get_color(numel(ind));
    v = [];
    for j = 1:numel(ind)-1 % assign 'speed' to each molecule
        dx0 = q*(xf_all(ind(j)) - xf_all(ind(j+1))).^2;
        dy0 = q*(yf_all(ind(j)) - yf_all(ind(j+1))).^2;
        dz0 = q*(zf_all(ind(j)) - zf_all(ind(j+1))).^2;
        v(j) = (dx0+dy0+dz0).^0.5/ (xpo);
    end
    trajec(i).v = mean(v);
    if abs(mean(v) - mean(mmax)) <= diff(mmax)
        plot3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind));
    end
    vs(i) = mean(v);
%     dx = [dx;q*(xf_all(ind)-mean(xf_all(ind)))];
%     dy = [dy;q*(yf_all(ind)-mean(yf_all(ind)))];
%     dz = [dz;q*(zf_all(ind)-mean(zf_all(ind)))];
%     ddz = [ddz; q*((diff(zf_all(ind))))./(diff(framenumber(ind))*0.02)];
%     mz = [mz; mean(diff(zf_all(ind)))*q];
%     mf = [mf; mean((framenumber(ind)))];
%     zd = zd + sum(zf_all(ind));
%     N{i} = fits(ind,3);
    n = [n;fits(ind,3) - mean(fits(ind,3))];
%     plot3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind),'LineWidth',1);
%     si = scatter3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind),-scl*llv(ind)./fits(ind,3),rgb);
%     si.MarkerFaceColor = si.MarkerEdgeColor;
%     plot3(q*(xf_all(ind)-xf_all(ind(1))),q*(yf_all(ind)-yf_all(ind(1))),i+q*(zf_all(ind)-zf_all(ind(1))))
    hold on
end
hold off
drawnow
figure
histogram(vs)
hold on
plot([min(mmax),min(mmax)],[0,100],'r');
plot([max(mmax),max(mmax)],[0,100],'r');
hold off
% clear uf uz
% uf = unique(mf);
% for i = 1:numel(uf)
%     ind = mf == uf(i);
%     uz(i) = mean(mz(ind));
% end
% 
% zdrift = gausssmooth(uz,145,180);
% sz = spline(uf,zdrift,framenumber);
% zfs = zf_all*q - sz;
% close all
% s = scatter3(xf_fixed*q,yf_fixed*q,zf_all*q,10,framenumber);
% s.MarkerFaceColor = s.MarkerEdgeColor;
% colormap('jet')
% axis equal
% 
% drawnow
% dz = [];
% dx = [];
% dy = [];
% hold on
% for i = 1:numel(trajec)
%     ind = trajec(i).t;
%     nums(i) = numel(ind);
%     indy = snr(ind) < msnr;
%     ind(indy) = [];
%     rgb = get_color(numel(ind));
%     dx = [dx;q*(xf_all(ind)-mean(xf_all(ind)))];
%     dy = [dy;q*(yf_all(ind)-mean(yf_all(ind)))];
%     dz = [dz;(zfs(ind)-mean(zfs(ind)))];
%     ddz = [ddz; q*((diff(zf_all(ind))))./(diff(framenumber(ind))*0.02)];
%     mz = [mz; mean(diff(zf_all(ind)))*q];
%     mf = [mf; mean((framenumber(ind)))];
%     zd = zd + sum(zf_all(ind));
%     N{i} = fits(ind,3);
%     n = [n;fits(ind,3) - mean(fits(ind,3))];
%     plot3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind),'LineWidth',1);
%     si = scatter3(q*xf_all(ind),q*yf_all(ind),q*zf_all(ind),-scl*llv(ind)./fits(ind,3),rgb);
%     si.MarkerFaceColor = si.MarkerEdgeColor;
%     plot3(q*(xf_all(ind)-xf_all(ind(1))),q*(yf_all(ind)-yf_all(ind(1))),i+q*(zf_all(ind)-zf_all(ind(1))))
%     hold on
% end
% xlim([10.4671 13.688])
% ylim([11.6372 14.9924])
% hold off
% axis equal
% vals(lk) = mean(abs(dz));
% end
% for i = 1:numel(N)
%     plot(N{i})
%     hold on
% end
% figure
% histogram(dz)
% fit_hist_gauss(dz)
