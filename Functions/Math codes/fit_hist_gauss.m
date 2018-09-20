function [fits, CovB] = fit_hist_gauss(X)
% this is a script to fit the histogram of N-D variable X w/ a gaussian and
% return the fitted gaussian values
h = histogram(X(:));
ys = h.Values;
xs = ((1:numel(ys)) - numel(ys)/2)*h.BinWidth;

% build the beta vector

beta0(1) = sum(ys.*xs)/sum(ys);
beta0(2) = max(ys);
beta0(3) = std(X(:));
beta0(4) = 0;
[fits, R, J,CovB] = nlinfit(xs,ys,@gauss1d,beta0);