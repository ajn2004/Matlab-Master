function yfit = gauss1d(beta,x)

yfit = (beta(2)*exp(- (x - beta(1)).^2./(2*beta(3)^2)) + beta(4)).';
end