function [y, t] = filter_cut(x,pad)
% Wrapper function to call filtcut gpu algorithm
y = [];
[m,n,o] = size(x);
tic
y = filtcut(x,pad,o);
t = toc;
end