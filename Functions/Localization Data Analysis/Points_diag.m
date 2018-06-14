close all;
clc;
clear Points x y
plot(xf_all,yf_all,'.');
% num = input('How Many Points are you selecting?');
[x,y] = ginput(2);
Points = {};
num = 1;

    
% num = 1;
for i = 1:num
    ind = xf_all >= min(x) & xf_all <= max(x);
    ind = ind & yf_all >= min(y) & yf_all <= max(y);
    ind = ind & abs(zf_all*133) <= 600;
    ind = ind & N <= 1000 & N >= 100;
    Points{i} = [xf_all(ind),yf_all(ind),zf_all(ind),framenum_all(ind),N(ind), xf_crlb(ind), yf_crlb(ind), zf_crlb(ind)];
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

 plot3(xf,yf,zf,'.')


