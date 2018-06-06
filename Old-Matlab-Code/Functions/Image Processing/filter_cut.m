function [y, t] = filter_cut(x,pad)
% Wrapper function to call filtcut gpu algorithm
y = [];
[m,n,o] = size(x);
type = class(x);
x = single(x);
tic
y = filtcut(x,pad,o);
t = toc;
if strcmp(type,'single')
    y = single(y);
end
end