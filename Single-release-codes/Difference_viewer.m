clearvars;
% close all;
clc;

% i1 = readtiff();
files = dir('local*');

for i = 1:numel(files)
    i1 = readtiff(files(i).name);
    i2(:,:,i) = i1(:,:,2)-i1(:,:,1);
i2(:,:,i) = (i2(:,:,i)>0).*i2(:,:,i);
imagesc(i2(:,:,i));
drawnow
waitforbuttonpress
end
i2 = (i2>0).*i2;
% imagesc(i2);