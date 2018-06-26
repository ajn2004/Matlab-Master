clearvars;
% close all;
clc;

% i1 = readtiff();
files = dir('*tif');

for i = 1:numel(files)
    i1 = readtiff(files(i).name);
    i2(:,:,i) = i1(:,:,2)-i1(:,:,1);
i2(:,:,i) = (i2(:,:,i)>0).*i2(:,:,i);
i3(:,:,i) = rollingball(i2(:,:,i));
% imagesc(i3(:,:,i));
% drawnow
% waitforbuttonpress
end
i2 = (i2>0).*i2;
% imagesc(i2);