function [z,params] = getdz(sigx,sigy,z_cal)
z =[];
x = (-0.5:0.001:0.5).';
z_cal = double(z_cal);
xs = z_cal(1);
gx = z_cal(2);
dx = z_cal(3);
ax = z_cal(4);
bx = z_cal(5);

ys = z_cal(6);
gy = z_cal(7);
dy = z_cal(8);
ay = z_cal(9);
by = z_cal(10);

% g = (gx+gy)/2;
% d = (dx+dy)/2;

sx = xs*(1 + ((x-gx)./dx).^2 + ax*((x-gx)./dx).^3 + bx*((x-gx)./dx).^4).^0.5; %defocusing stuff
sy = ys*(1 + ((x-gy)./dy).^2 + ay*((x-gy)./dy).^3 + by*((x-gy)./dy).^4).^0.5;

% sx = xs*(1 + ((x-g)./d).^2 + ax*((x-g)./d).^3 + bx*((x-g)./d).^4).^0.5; %defocusing stuff
% sy = ys*(1 + ((x+g)./d).^2 + ay*((x+g)./d).^3 + by*((x+g)./d).^4).^0.5;


for i = 1:numel(sigx)
    D = ((sigx(i).^0.5-sx.^0.5).^2 + (sigy(i).^0.5-sy.^0.5).^2).^0.5;
%     D = ((sigx(i)-sx).^2 + (sigy(i)-sy).^2).^0.5;
%     D = ((sigx(i).^0.5-sx.^0.5).^2 + (sigy(i).^0.5-sy.^0.5).^2);
    ind = find(D == min(D), 1);
    try
    z(i,1) = x(ind(1));  % emperically found to be the axial zoom factor by fitting a line between board movements and position measurements
    catch
        z(i,1) = -300000;
    end
end
params = [x,sx,sy];