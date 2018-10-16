% close all;
clc;
clear Points x y
% xf_a = xf_all - drifts(floor(framenum_all/3),1);
% yf_a = yf_all - drifts(floor(framenum_all/3),2);
xf_a = xf_all;
yf_a = yf_all;
plot(xf_a,yf_a,'.');
% num = input('How Many Points are you selecting?');
[x,y] = ginput(2);
Points = {};
num = 1;

    
% num = 1;
for i = 1:num
    ind = xf_a >= min(x) & xf_a <= max(x);
    ind = ind & yf_all >= min(y) & yf_all <= max(y);
    ind = ind & abs(zf_all*133) <= 600;
    ind = ind & N <= 500 & N >= 75;
    Points{i} = [xf_a(ind),yf_a(ind),zf_all(ind),framenum_all(ind),N(ind), xf_crlb(ind), yf_crlb(ind), zf_crlb(ind),cents(ind,:)];
end



% save('Points_frm_bump50.mat','Points');
data = Points{1};
Nf = data(:,5);
x = (data(:,4));
xf = data(:,1)*0.133*1000;
yf = data(:,2)*0.133*1000;
zf = data(:,3)*0.133*1000;
xf_e = data(:,6).^0.5*133;
yf_e = data(:,7).^0.5*133;
zf_e = data(:,8).^0.5*133;
cx = data(:,9);
cy = data(:,10);
 plot3(xf,yf,zf,'.')

dfc = ((xf-mean(xf)).^2 + (yf - mean(yf)).^2 + (zf - mean(zf)).^2).^0.5;
for i = 1 : numel(xf) - 1
    d3(i) = ((yf(i) - yf(i+1))^2+(xf(i) - xf(i+1))^2+(zf(i) - zf(i+1))^2)^0.5;
    d2(i) = ((yf(i) - yf(i+1))^2+(xf(i) - xf(i+1))^2)^0.5;
end
hold on
