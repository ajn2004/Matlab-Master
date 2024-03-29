% Get the signal to noise
close all;
clearvars
i1 = readtiff()/33.33; % 33.33 is the measured gain on my camera
pix = 3;
[m,n,o] = size(i1);
for i = 1:o-1
i2(:,:,i) = i1(:,:,i+1) - i1(:,:,i);
snr1(:,:,i) = (i2(:,:,i).*(i2(:,:,i)>0))./i1(:,:,i).^0.5;
end

imagesc(i2(:,:,1))
[x,y] = ginput(1);
x = round(x);
y = round(y);
wind = -pix:pix;
sub1 = i1(y+wind, x+wind,1);
sub2 = i2(y+wind, x+wind,1);
snr = sub2./sub1.^0.5;
figure
imagesc(snr);
% snr1 = (i2(:,:,1).*(i2(:,:,1)>0))./i1(:,:,1).^0.5;
figure
% imagesc(snr1);
set_scale(mean(snr(:,:,1:98),3), 0.133, 4);