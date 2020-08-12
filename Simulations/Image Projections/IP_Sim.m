%% Imaging Projection results
% Create PSFs w/ varying behavior, use various projections to show their
% impact on different imaging conditions
% AJN 5/12/19
clearvars; close all; clc;

sx = 1.3;
N = 3000;
O = 1.5;

z = -1:0.02:1;
dx = 0.3273;
s = sx*(1 + ((z)./dx).^2).^0.5;

%% Build constant low signal high noise molecules
for i = 1:100
    i1 = get_psf(0,0,sx,sx);
    B = ((N/10)^2-N)/numel(i1(:));
    i2(:,:,i) = double(imnoise(uint16(N*i1 + B),'poisson'));
    imagesc(i2(:,:,i),[200,1000]);
    axis image
    colormap('gray');
    M(i) = getframe(gcf);
end
writetiff(mean(i2,3),'low_snr_average.tif')
writetiff(max(i2,[],3),'low_snr_maximal.tif')
writetiff(std(i2,[],3),'low_snr_std.tif')
movie2gif(M,'lowsnr_molecule.gif','DelayTime',0.04,'LoopCount',Inf);

%% Build high signal to noise varying width molecules
for i = 1:100
    i1 = get_psf(0,0,s(i),s(i));
    B = 0;
    i2(:,:,i) = double(imnoise(uint16(N*i1 + B),'poisson'));
    imagesc(i2(:,:,i),[0,300]);
    axis image
    colormap('gray');
    M(i) = getframe(gcf);
end
writetiff(mean(i2,3),'blur_average.tif')
writetiff(max(i2,[],3),'blur_maximal.tif')
writetiff(std(i2,[],3),'blur_std.tif')
movie2gif(M,'lowsnr_molecule.gif','DelayTime',0.04,'LoopCount',Inf);

