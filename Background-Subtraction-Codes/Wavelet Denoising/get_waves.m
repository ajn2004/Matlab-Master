function [W, i3] = get_waves(i1,k)
% This function will get the wavelet transform W and corresponding smoothed
% images I of i1 for k levels
[m,n] = size(i1);

lvls = k;
baselet = [1/16, 1/4, 3/8, 1/4, 1/16];

i3(:,:,1) = double(i1);
for i = 1:lvls % loop over K levels
    wvlt = [baselet(1), zeros(1,2^(i-1)-1),baselet(2), zeros(1,2^(i-1)-1),baselet(3), zeros(1,2^(i-1)-1),baselet(4), zeros(1,2^(i-1)-1),baselet(5)]; % construct wavelet
    % perform seperable convolution
    for j = 1:m % convolve all rows
        i3(j,:,i+1) = conv(i3(j,:,i),wvlt,'same');
    end
    for j = 1:n % convolve all columns
        i3(:,j,i+1) = conv(i3(:,j,i+1),wvlt,'same');
    end
%     hold on
%     plot(wvlt)
    W(:,:,i) = i3(:,:,i) - i3(:,:,i+1); % wavelet plane is defined as difference of smoothed images
end


