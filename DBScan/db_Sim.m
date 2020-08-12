% DB Explainer
close all
points = 3000;
maxx = 4;
sfold = 'C:\Users\AJN Lab\Dropbox\Lab Meetings AJN\July 22 2019\';

xr = rand(points/10,1)*4;
yr = rand(points/10,1)*4;

cen = 2;
r = 1;

theta = rand(points,1)*2*pi;

xc = r*cos(theta) + rand(points,1)*0.2 + cen;
yc = r*sin(theta) + rand(points,1)*0.2 + cen;

f = plot(xc,yc,'.k')
hold on
ax = gca;
plot(ax,xr,yr,'.k')
title('Raw Data')
saveas(f,[sfold,'raw_db_sim.tif']);
xf = [xc;xr];
yf = [yc;yr];

clust = DB_scan([xf,yf],0.1,4);
id = clust == 1;
plot(xf,yf,'.')
hold off
g = plot(xf(id),yf(id),'.')
ax = gca;
hold on
plot(xf(logical(1-id)),yf(logical(1-id)),'.k');
legend('Cluster','Not Cluster')
title('Cluster Assignment')
saveas(g,[sfold,'clustered_DB_SIM.tif']);