% clc; close all; clearvars;
pix = 100;
mols = 20;
bness = 200;
bkn = sqrt(bness);
i1 = zeros(pix);

for i = 1:poissrnd(mols)
    l = randi(pix);
    k = randi(pix);
    i1(l,k) = poissrnd(bness);
end
[x,y] = meshgrid(-5:5,-5:5);
z = exp(-(x.^2+y.^2)/(2*1.6^2));
z = z/sum(z(:));
i2 = gpu_conv(i1,z);
imagesc(i2);
i3 = imnoise(uint16(i2 + bkn),'Poisson');
imagesc(i3)
i4 = bandpass(i3);
i5 = rollingball(i3);
subplot(2,2,1);surf(i4); axis image; title('BKNG sub');
% subplot(2,2,2);imagesc(i3); axis image; title('Noised');
% subplot(2,2,3);imagesc(i2); axis image; title('Simmed');

ii = i2(6:end-6,6:end-6)-i4(6:end-6,6:end-6);
ij = i2(6:end-6,6:end-6)-i5(6:end-6,6:end-6);
ik = i4(6:end-6,6:end-6)-i5(6:end-6,6:end-6);
subplot(2,2,2);histogram(ik(:));  title('Compare');
subplot(2,2,3);histogram(ij(:)); title('rollingball');
subplot(2,2,4);histogram(ii(:));  title('bandpass');
