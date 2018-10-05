% LP filter comparison
va = [0,30];

i1 = readtiff();
im = 200;

i2 = i1(:,:,im);
[m,n,] = size(i2);

for i = 1:max([m,n])
    i3 = lp_filt(i2,i*2);
    subplot(1,2,1);
    imagesc(i3);
    axis image
    subplot(1,2,2);
    surf(i3);
    zlim([min(i2(:)),max(i2(:))]);
    view(va);
    M(i) = getframe(gcf);
end