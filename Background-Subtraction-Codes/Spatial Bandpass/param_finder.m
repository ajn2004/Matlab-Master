% parameter testing
% this code will test the parameter
% close all; clc;
% clearvars;
im1 = readtiff('cell4 488.tif');
i1 = im1(:,:,1);
% for i = 1:100
%     for j = 1:100
        i2 = bandpass(i1);
%         i3(:,:,i,j) = i2;
imagesc(i2);
axis image
drawnow
%         M((i-1)*100 + j) = getframe(gcf);
%     end
% end