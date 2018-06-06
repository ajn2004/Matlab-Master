% PSF_timeline
% This will allow a use to select a PSF from a summed frame and cut out the
% corresponding region on all frames, then save that as a separate tiff
% function psf_cut(num)
psfs = 1;
num = 2;
[fname, fpath] = uigetfile('*tif');
pixw = 10;
im1 = readtiff(fname);
sim1 = sum(im1,3);
imagesc(sim1, [min(sim1(:)), 0.75*max(sim1(:))])
colormap('gray')
axis image
[x,y] = ginput(psfs);
% Center around the maximal pixel of the sum
[row, col] = find(sim1 == max(max(sim1(round(y)-pixw:round(y)+pixw,round(x)-pixw:round(x)+pixw))));
[m,n,o] = size(im1);
im2 = rollingball(im1);
% imname = 'psf_';
A = cell(psfs,1);
% A=zeros(pixw,pixw,o);
for i = 1:numel(x)
    A{i} = im2(round(y(i))-pixw:round(y(i))+pixw,round(x(i))-pixw:round(x(i))+pixw,:);
%     writetiff(A{i},[imname,num2str(num)]);
end
save(['PSF_',num2str(num)],'A');