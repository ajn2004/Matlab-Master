function z_cal = get_z_params(x, sigx, sigy)
% This function takes in z space, sigx_all, and sigy_all and fits these
% values to the defocusing function returning the defocusing parameters
% this will return z_cal in the format of [sx, gx, dx, ax, bx, sy, gy, dy, ay, by]
x = double(x);
sigx = double(sigx);
sigy = double(sigy);
fun = @(beta1,x)beta1(1)*(1+((x-beta1(2))./beta1(3)).^2 + beta1(4)*((x-beta1(2))./beta1(3)).^3 + beta1(5)*((x-beta1(2))./beta1(3)).^4).^0.5;

x0 = [1.2, 0, 0.5, 0, 0];
fx = lsqcurvefit(fun,x0,x(:),sigx(:));
fy = lsqcurvefit(fun,x0,x(:),sigy(:));

z_cal = [fx,fy];

% yx = z_cal_fit(x,z_cal(1:5));
% yy = z_cal_fit(x,z_cal(6:end));
% close all
% plot(x,sigx,'.');
% hold on
% plot(x,sigy,'.');
% plot(x,yx)
% plot(x,yy)
% drawnow
