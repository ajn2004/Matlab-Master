clearvars;
close all;
clc;

load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat');
% [z,params] = getdz(1,1,cal.red.z_cal);
% f = figure('Units','Normalized','Outerposition',[0 0 1 1])
[x,y] = meshgrid(-8:8);
zs = cal.red.z0s;
sx = cal.red.sx;
sy = cal.red.sy;
ind = abs(zs) <= 0.8;

zs1 = zs(ind);
sx1 = sx(ind);
sy1 = sy(ind);
zs = -0.8:0.01:0.8;
sx = interp1(zs1,sx1,zs);
sy = interp1(zs1,sy1,zs);
sx = sx(1:end-1);
zs = zs(1:end-1);
sy = sy(1:end-1);
% sx = 3*zs.^2 +1.2 ;
% sy = sx;
count = 1;
for i = 1:1:numel(sx)
    subplot(3,6,[1,7]);
    plot(0,zs(i),'.k','markersize',10)
    ylim([min(zs),max(zs)])
    ylabel('Axial Position of Molecule')
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    z = exp(-x.^2/(2*sx(i).^2) - y.^2/(2*sy(i).^2));
    subplot(3,6,[2:6,8:12]);
    imagesc(z)
    set(gca,'visible','off')
    hold on
    plot(9+[-2*sx(i),2*sx(i)],[9,9],'r')
    plot(9+[-2*sx(i),-2*sx(i)],[8.5,9.5],'r')
    plot(9+[2*sx(i),2*sx(i)],[8.5,9.5],'r')
    plot([9,9],9+[-2*sy(i),2*sy(i)],'g')
    plot([8.5,9.5],9-[2*sy(i),2*sy(i)],'g')
    plot([8.5,9.5],9+[2*sy(i),2*sy(i)],'g')
    
    hold off
    title('Camera Image')
    xlabel('Position in Pixels')
    drawnow
    axis image
    subplot(3,6,13:18)
    plot(zs,2*sx,'r','linewidth', 2)
    hold on
    plot(zs,2*sy,'g','linewidth', 1)
    plot([zs(i),zs(i)],[0,2*max([max(sx),max(sy)])],'k')
    xlim([-0.8, 0.8])
    ylim(2*[1 4])
    hold off
    xlabel('Axial Position')
    ylabel('Gaussian Width')
    legend('2\sigma_X','2\sigma_Y','location','best')
    title('Widths versus Axial Position')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    M(count) = getframe(gcf);
    count = count +1 ;
end

for i = numel(sx):-1:1
    subplot(3,6,[1,7]);
    plot(0,zs(i),'.k','markersize',10)
    ylim([min(zs),max(zs)])
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    ylabel('Axial Position of Molecule')
    z = exp(-x.^2/(2*sx(i).^2) - y.^2/(2*sy(i).^2));
    subplot(3,6,[2:6,8:12]);
    imagesc(z)
    set(gca,'visible','off')
    hold on
    plot(9+[-2*sx(i),2*sx(i)],[9,9],'r')
    plot(9+[-2*sx(i),-2*sx(i)],[8.5,9.5],'r')
    plot(9+[2*sx(i),2*sx(i)],[8.5,9.5],'r')
    plot([9,9],9+[-2*sy(i),2*sy(i)],'g')
    plot([8.5,9.5],9-[2*sy(i),2*sy(i)],'g')
    plot([8.5,9.5],9+[2*sy(i),2*sy(i)],'g')
    
    hold off
    title('Camera Image')
    xlabel('Position in Pixels')
    drawnow
    axis image
    subplot(3,6,13:18)
    plot(zs,2*sx,'r','linewidth', 2)
    hold on
    plot(zs,2*sy,'g','linewidth',1)
    plot([zs(i),zs(i)],[0,2*max([max(sx),max(sy)])],'k')
    xlim([-0.8, 0.8])
    ylim(2*[1 4])
    hold off
    xlabel('Axial Position')
    ylabel('Gaussian Width')
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    legend('2\sigma_X','2\sigma_Y','location','best')
    title('Widths versus Axial Position')
    M(numel(M)+1) = getframe(gcf);
end

degs = 5;
sx = sx - mean(sx);
sy = sy - mean(sy);
[X Y Z] = rot_mat(deg2rad(degs));
data = [sx;sy];
count = 1;
% for i = 0:degs:360
    
%     data = Z(1:2,1:2)*data;
    plot3(data(1,:),data(2,:),zs,'linewidth',4)
    xlim([min(sx)*1.2,max(sx)*1.2])
    ylim([min(sx)*1.2,max(sx)*1.2])
    zlim([min(zs)*1.1, max(zs)*1.1])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    xlabel('Guassian Width - X')
    ylabel('Gaussian Width - Y')
    zlabel('Axial Position')
%     M(count) = getframe(gcf);
%     count = count +1;
% end