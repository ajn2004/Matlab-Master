%% Side-by-side movies
% Play two tiff files side by side and record the frames for a subsequent
% figure
clearvars; close all; clc;

% Select images and load them into memory
i1 = readtiff();
i2 = readtiff();

% get size of images
[m1,n1,o1] = size(i1);
[m2,n2,o2] = size(i2);

%% Make frames
% loop over the minimum number of frames
for i = 1:min([o1,o2])
    
    % Display first image on the left
    subplot(1,2,1); 
    imagesc(i1(:,:,i), [min(i2(:))*1.3, 0.8*max(i2(:))]);
    title('Image 1');
    axis image
    
    % Display second image on the right
    subplot(1,2,2);
    imagesc(i2(:,:,i), [min(i2(:))*1.3, 0.8*max(i2(:))]);
    axis image
    title('Image2');
    M(i) = getframe(gcf);
end

