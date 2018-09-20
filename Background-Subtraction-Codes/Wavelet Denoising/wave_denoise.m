% Wavelet denoising
% This script is being done in the spirit of Smal et al 2010 IEEE
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5109713
% AJN 9/18/18
clearvars; close all
lvls = 4; % number  of K wavelet images to create
try % try loading an image
    i1 = readtiff('local_2.tif');
catch
    i1 = readtiff();
end
% i1 is a loaded image

ind = 2;
% i2 = i1(:,:,ind)-i1(:,:,1); % grab the image w/ a molecule
i2 = double(i1(:,:,ind));
% Figure tabset
f = figure;
tg = uitabgroup(f);
t1 = uitab(tg,'title','Original Image');
ax = axes(t1);
surf(ax,i2)
% axis image

baselet = [1/16, 1/4, 3/8, 1/4, 1/16];
% plot(baselet)
[m,n] = size(i2);
i3 = zeros(m,n,lvls+1);
% W = i3;
% for l = 2
%     i2 = i1(:,:,l);
i3(:,:,1) = i2;
% W(:,:,1) = i2;
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
t2 = uitab(tg,'title','Reconstructed Image');
ax = axes(t2);
imagesc(ax,sum(W,3) + i3(:,:,end) - i2)% This is the inverse transform represented as an image, it should be identical to i2
% title('reconstructed image')
% axis image

% Denoise the system via thresholding
for i = 2:lvls+1
    ii3 = W(:,:,i);
    thrsh = std(ii3(:));   
    im3 = W(:,:,i).^2-3*thrsh^2;
    i4(:,:,i) = (W(:,:,i).^-1).*(im3.*(im3>0));
%     i3(:,:,i) = i3(:,:,i).*i3(:,:,i)>thrsh;
end
t3 = uitab(tg,'title','Denoised Image');
ax = axes(t3);
surf(ax,sum(W(:,:,2:i),3) + i4(:,:,i))% This is the inverse transform represented as an image, it should be identical to i2
% axis image
den = sum(W(:,:,2:i),3) + i4(:,:,i);
i5(:,:,l) = den;
% end
im1 = i1(:,:,2) - den;
figure
i2 = (i5(:,:,2)-i5(:,:,1));
