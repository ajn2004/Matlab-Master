function [fluor, xs] = func_measure_regions(im1, cents, psize)
% take in image, centers, and region size, return time average

%preallocate
fluor = zeros(1,numel(im1(1,1,:)));
% im1s = zeros(psize*2+1,psize*2+1,numel(im1(1,1,:)),numel(cents(:,1)));
% average over square regions of size -psize:psize
for i = 1:numel(cents(:,1))
    im1s(:,:,:,i) = im1(cents(i,2)-psize:cents(i,2)+psize, cents(i,1)-psize:cents(i,1)+psize, :);
end
count =1;
for i = 1:numel(im1s(1,1,1,:))
    xs(count,:) = mean(mean(im1s(:,:,:,i),1),2);
    count = count +1;
end
fluor = mean(mean(mean(im1s,1),2),4);