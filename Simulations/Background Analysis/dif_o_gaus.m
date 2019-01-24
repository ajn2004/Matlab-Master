% Testing a difference in gaussians filters
close all
clearvars
clc
pixw = 4;
scle = 4;
ywind = 190 + (-(scle*pixw):(scle*pixw));
xwind = 233 + (-(scle*pixw):(scle*pixw));
sig = 1;
for i = 1:10
relscale = 1+i/10;
[xp,yp] = meshgrid(-pixw:pixw);
i1 = readtiff('C:\Users\AJN Lab\Dropbox\Data\10-30-18 vglut-ph\local2\local_27.tif');
g1 = exp(-(xp.^2+yp.^2)/(2*sig^2));
g1 = g1/sum(g1(:));
g2 = exp(-(xp.^2+yp.^2)/(2*(relscale*sig)^2));
g2 = g2/sum(g2(:));
im2 = i1(:,:,2);
im1 = i1(:,:,1);
figure('Units','Normalized','Outerposition',[0 0.5 0.5 0.5])
dim = im2-im1;
imagesc(dim(ywind,xwind));
% axis equal
title('Original')
figure('Units','Normalized','Outerposition',[0.5 0.5 0.5 0.5])
dim2 = conv2(im2,g1,'same')-conv2(im1,g2,'same');
imagesc(dim2(ywind,xwind));
title('Diff of convs')
figure('Units','Normalized','Outerposition',[0 0 0.5 0.5])
dim3 = conv2(im2-im1,g1-g2,'same');
imagesc(dim3(ywind,xwind))
title('g1-g2')
figure('Units','Normalized','Outerposition',[0.5 0 0.5 0.5])
dim4 = conv2(im2,g1,'same')-conv2(im1,g2,'same') - conv2(im2-im1,g1-g2,'same');
imagesc(dim4(ywind,xwind));
title('Diff of methods')
pause(1)
% close all
end