close all;
clc;
clear Points x y
plot(xf_all,yf_all,'.');
% num = input('How Many Points are you selecting?');
[x,y] = ginput(2);
Points = {};
num = 1;
if x(1) < x(2)
    xx = x(2);
    x(2) = x(1);
    x(1) = xx;
    clear xx
end
if y(1) < y(2)
    yy = y(2);
    y(2) = y(1);
    y(1) = yy;
end
    
% num = 1;
for i = 1:num
    ind = N >=20& xf_all > x(2) & xf_all < x(1) & yf_all < y(1) & yf_all > y(2);
    Points{i} = [xf_all(ind),yf_all(ind),zf_all(ind),framenum_all(ind),N(ind), xf_crlb(ind), yf_crlb(ind), zf_crlb(ind)];
    if sum(ind) > 5
        plot(framenum_all(ind),zf_all(ind));
    end
    hold on
end
hold off
legend('1','2','3','4');


% clearvars -except Points
for i = 1:numel(Points)
p = i;
% save('Points_frm_bump50.mat','Points');
data = Points{p};
% sigx = data(:,6);
% sigy = data(:,7);
Nf = data(:,5);
x = (data(:,4));
xf = data(:,1)*0.133*1000;
yf = data(:,2)*0.133*1000;
zf = data(:,3)*0.133*1000;
xf_e = data(:,6).^0.5*133;
yf_e = data(:,7).^0.5*133;
zf_e = data(:,8).^0.5*133;
plot(xf,zf,'.');
[xx,z] = ginput(2);
if z(1) > z(2)
    zz = z(1);
    z(2) = z(1);
    z(1) = zz;
end
% ind = zf <= z(2) & zf >= z(1);
% zf = zf(ind);
% x = x(ind);
% plot(x,zf)
% ind = x > 2000 & x < 3000;
% % plot(x(ind),zf(ind))
% plot(x*0.01,zf);
% % s = scatter3(xf,yf,zf,[],x);
% xlabel('Time [s]');
% ylabel('Axial Position[nm]')
hold on
end
hold off
% plot(x,0.2413*x-874.2072,'r');
% plot(x,2.4594*10^-5*x.^2+0.0707*x - 594.9947,'g')
% hold off
% xbins = 0.008*[3619 , 3797; 3808, 4002; 4005, 4202; 4223, 4381; 4397, 4566; 4587, 4803; 4811,5013;  5048, 5219; 5245, 5405; 5429, 5603; 5632, 5797; 5813, 5968; 5991, 6198; 6215,6422; 6443, 6598; 6626, 6792; 6834,6994; 7018, 7210; 7238, 7384; 7396,7603; 7626, 7809; 7820, 8003;8043,8208];       
% xbins =[0, 2046; 2100, 2245; 2276, 2388; 2424, 2610; 2645, 2799; 2821, 2998; 3035, 3129; 3142, 3382; 3399, 3556; 3580, 3710; 3760, 3881; 3906,4055;4068,4191;4210,4355;4368,4499;4510,4615;4638,4757;4781,4888];
% for i = 1:numel(xbins(:,1))
%     ind = x > xbins(i,1) & x < xbins(i,2);
%     mx(i) = mean(x(ind));
%     mea(i) = mean(zf(ind));
%     stdd(i) = std(zf(ind));
%     val(i) = sum(ind);
%     clear ind
% end
% errorbar(mx,mea,stdd./val.^0.5,'r.')
% a = polyfit(mx,mea,1);
% plot(mx,a(1)*mx+a(2),'r')
% hold off
% plot([4219, 4219],[min(zf), max(zf)],'r');
% plot([4355, 4355],[min(zf), max(zf)],'r');
% plot([3142, 3142],[min(zf), max(zf)],'r');
% plot([3359, 3359],[min(zf), max(zf)],'r');
% hold off
% xlabel('Frame Number');
ylabel('Fitted Position[nm]');
title('Bump Test of 50 nm');
% figure
% histogram(N)
% x = (data(:,4)-230)*0.005;
% z_cal = get_z_params(x,sigx,sigy);
% xlabel('Axial Position [nm]');
% ylabel('Sigma width [pixels]');
% title('Axial Calibration and Fit');
