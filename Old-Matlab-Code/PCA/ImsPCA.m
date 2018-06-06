function [yout, t] = ImsPCA(i1,pad,varargin)
if numel(varargin) == 1
    k = varargin{1};
elseif numel(varargin) ==0
    k = 1;
end
[m,n,o] = size(i1);
yout = i1;
for i = 1:k:o
    try
    [ys, t(i:i+k-1)]  = IPCA(i1(:,:,i:i+k-1),pad);
    catch
        [ys, t(i:i+k-1)]  = IPCA(i1(:,:,i:end),pad);
    end
    yout(:,:,i:i+numel(ys(1,1,:))-1) = ys;
end
end