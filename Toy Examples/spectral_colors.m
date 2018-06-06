%% Spectral Colors
% a code that creates a spectral color bar based on gaussian distributions over
% R G B space, builds a resulting movie as these weights shift

spc=10;
x = 0:0.01:1;

for j = 1:100
    r = exp(-(x-0).^2/(2*.25^2));
    g = exp(-(x-(j-1)/100).^2/(2*.25^2));
    b = exp(-(x-1).^2/(2*.25^2));
    rgb = [r;g;b].';
    im1 = zeros(20,spc*numel(rgb(:,1)),3);
    for i = 1:numel(rgb(:,1))
        im1(:,(i-1)*spc+1:i*spc,1) = rgb(i,1);
        im1(:,(i-1)*spc+1:i*spc,2) = rgb(i,2);
        im1(:,(i-1)*spc+1:i*spc,3) = rgb(i,3);
    end
    subplot(3,4,[1,2,3,5,6,7,9,10,11]);imagesc(im1)
    subplot(3,4,4);plot(x,r,'r'); title('Red Value');
    subplot(3,4,8);plot(x,g,'r'); title('Green Value');
    subplot(3,4,12);plot(x,b,'r'); title('Blue Value');
    drawnow
    M(j) = getframe(gcf);
end
for k = 1:100
    r = exp(-(x-0).^2/(2*.25^2));
    g = exp(-(x-(100-k)/100).^2/(2*.25^2));
    b = exp(-(x-1).^2/(2*.25^2));
    rgb = [r;g;b].';
    im1 = zeros(20,spc*numel(rgb(:,1)),3);
    for i = 1:numel(rgb(:,1))
        im1(:,(i-1)*spc+1:i*spc,1) = rgb(i,1);
        im1(:,(i-1)*spc+1:i*spc,2) = rgb(i,2);
        im1(:,(i-1)*spc+1:i*spc,3) = rgb(i,3);
    end
    subplot(3,4,[1,2,3,5,6,7,9,10,11]);imagesc(im1)
    subplot(3,4,4);plot(x,r,'r'); title('Red Value');
    subplot(3,4,8);plot(x,g,'r'); title('Green Value');
    subplot(3,4,12);plot(x,b,'r'); title('Blue Value');
    drawnow
    M(j+k) = getframe(gcf);
end