function [fits, CovB] = fit_hist_gauss(X)
% this is a script to fit the histogram of N-D variable X w/ a gaussian and
% return the fitted gaussian values
[N,edge] = histcounts(X(:));
ys = N;
xs = (edge(1:end-1))+(edge(2)-edge(1))/2;
bar(xs,ys)
% build the beta vector

beta0(1) = sum(ys.*xs)/sum(ys);
beta0(2) = max(ys);
beta0(3) = std(X(:));
beta0(4) = 0;
[fits, R, J,CovB] = nlinfit(xs,ys,@gauss1d,beta0);
xt = min(xs):0.001:max(xs);
yfit = (fits(2)*exp(- (xt - fits(1)).^2./(2*fits(3)^2)) + fits(4));
hold on
plot(xt,yfit,'r')
hold off