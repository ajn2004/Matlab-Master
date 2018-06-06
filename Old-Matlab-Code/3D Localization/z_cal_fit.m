function [y] = z_cal_fit(x,beta1)

y = beta1(1)*(1+((x-beta1(2))./beta1(3)).^2 + beta1(4)*((x-beta1(2))./beta1(3)).^3 + beta1(5)*((x-beta1(2))./beta1(3)).^4).^0.5;
end