function [yout, t] = ImsPCA(i1,pad,varargin)
if numel(varargin) == 1
    k = varargin{1};
elseif numel(varargin) ==0
    k = 1;
end
[m,n,o] = size(i1);
yout = i1;
count = 1;
for i = 1:k:o
    try
    tic
        [ys, t(i:i+k-1)]  = IPCA(i1(:,:,i:i+k-1),pad);
    t(count) = toc;
    count = count +1;
%     ajn_wait(mean(t), i,o/k);
%     disp('Calculating PCA');
    catch
        [ys, t(i:i+k-1)]  = IPCA(i1(:,:,i:end),pad);
    end
    yout(:,:,i:i+numel(ys(1,1,:))-1) = ys;
end
end