function rgb = get_color(divs)
x = (1/divs:1/divs:1) - 1/divs;
spc = 2;
r = exp(-(x-0).^2/(2*.25^2));
g = exp(-(x-0.5).^2/(2*.25^2));
b = exp(-(x-1).^2/(2*.25^2));
rgb = [r;g;b].';
% figure
% im1 = zeros(5,spc*numel(rgb(:,1)),3);
%     for i = 1:numel(rgb(:,1))
%         im1(:,(i-1)*spc+1:i*spc,1) = rgb(i,1);
%         im1(:,(i-1)*spc+1:i*spc,2) = rgb(i,2);
%         im1(:,(i-1)*spc+1:i*spc,3) = rgb(i,3);
%     end
%     imagesc(im1)
end