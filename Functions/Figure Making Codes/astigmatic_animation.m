clearvars;
% close all;
clc;

load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat');
[z,params] = getdz(1,1,cal.z_cal);
% f = figure('Units','Normalized','Outerposition',[0 0 1 1])
[x,y] = meshgrid(-8:8);
zs = params(:,1);
sx = params(:,2);
sy = params(:,3);
count = 1;
for i = 1:18:numel(params(:,1))
    subplot(3,6,[1,7]);
    plot(0,zs(i),'.')
    ylim([min(zs),max(zs)])
    ylabel('Axial Position in um')
    z = exp(-x.^2/(2*sx(i).^2) - y.^2/(2*sy(i).^2));
    subplot(3,6,[2:6,8:12]);
    imagesc(z)
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
    plot(zs,2*sx)
    hold on
    plot(zs,2*sy)
    plot([zs(i),zs(i)],[0,2*max([max(sx),max(sy)])],'r')
    xlim([-0.8, 0.8])
    ylim(2*[1.1 3.5])
    hold off
    xlabel('Axial Position in um')
    ylabel('Width in pixels')
    legend('Sigma X','Sigma Y')
    title('Calibration Curve for Astigmatism')
    M(count) = getframe(gcf);
    count = count +1 ;
end

for i = numel(params(:,1)):-18:1
    subplot(3,6,[1,7]);
    plot(0,zs(i),'.')
    ylim([min(zs),max(zs)])
    ylabel('Axial Position in um')
    z = exp(-x.^2/(2*sx(i).^2) - y.^2/(2*sy(i).^2));
    subplot(3,6,[2:6,8:12]);
    imagesc(z)
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
    plot(zs,2*sx)
    hold on
    plot(zs,2*sy)
    plot([zs(i),zs(i)],[0,2*max([max(sx),max(sy)])],'r')
    xlim([-0.8, 0.8])
    ylim(2*[1.1 3.5])
    hold off
    xlabel('Axial Position in um')
    ylabel('Width in pixels')
    legend('Sigma X','Sigma Y')
    title('Calibration Curve for Astigmatism')
    M(numel(M)+1) = getframe(gcf);
end