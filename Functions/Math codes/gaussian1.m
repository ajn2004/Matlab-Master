function y = gaussian1(x,beta0)
y = beta0(1)*exp(-(x-beta0(2)).^2/(2*beta0(3)^2))+beta0(4);