function y = gaussian(sig,pixw)
% Return a normalized gaussion of size 2*pixw+1 w/ a std. deviation of sig
x = -pixw:pixw;

y = exp(-x.^2/(2*sig));
y = y/sum(y);