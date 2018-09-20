function [W] = thresh_waves(W)
[m,n,o] = size(W);
lvls = o;
% Thresholding the components
% Thresholding is done to the W matrix
for i = 1:lvls
    ii3 = W(:,:,i);
    thrsh = std(ii3(:));
    im3 = ii3.^2-9*thrsh^2;
    W(:,:,i) = (ii3.^-1).*(im3.*(im3>0));
end