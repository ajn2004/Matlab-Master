function [z] = getdz(sigx,sigy)
z =[];
x = (-1:0.001:1).';
% load('C:\Users\AJN Lab\Desktop\Code Development\Matlab Testing Folder\3D Localization\Z calibrations\3d_theta.mat');
load('C:\Users\AJN Lab\Desktop\Code Development\Matlab Testing Folder\3D Localization\Z calibrations\z_cal.mat');
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


sx = xs*(1 + ((x-gx)./dx).^2 + ax*((x-gx)./dx).^3 + bx*((x-gx)./dx).^4).^0.5; %defocusing stuff
sy = ys*(1 + ((x-gy)./dy).^2 + ay*((x-gy)./dy).^3 + by*((x-gy)./dy).^4).^0.5;
% ss = func_eval_NN(x,theta1,theta2);

% sx = ss(:,1)*ysc(1) + ymi(1);
% sy = ss(:,2)*ysc(2) + ymi(2);

for i = 1:numel(sigx)
    D = ((sigx(i).^0.5-sx.^0.5).^2 + (sigy(i).^0.5-sy.^0.5).^2);
    ind = find(D == min(D), 1);
    try
    z(i,1) = 1.0015*x(ind(1));  % emperically found to be the axial zoom factor by fitting a line between board movements and position measurements
    catch
        z(i,1) = -300000;
    end
end