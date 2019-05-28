function yfit = gauss1d(beta0,x)

yfit = (beta0(2)*exp(- (x - beta0(1)).^2./(2*beta0(3)^2)) + beta0(4)*0);
end