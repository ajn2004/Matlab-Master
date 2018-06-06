function [yout, t] = IPCA(i1, pad)
% Function that performs an image filter in the form of Principle Component
% Analysis

[m,n,o] = size(i1);

[y, t] = filter_cut(i1,pad);
[V, D] = eig(cov(y - mean(y)));
v = V(:,end);
y1 = v.'*(y- mean(y)).';
yout = reshape(y1,m,n,o);
end