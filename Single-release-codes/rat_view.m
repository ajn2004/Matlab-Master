function [thrsh] = rat_view(i1)
% this is a function to find molecules in noisey images by looking at the
% ratio of a 'signal' image and the 'background' image
% close all
[m,n,o] = size(i1);
maxv = reshape(max(max(i1)),1,o);
minv = reshape(min(min(i1)),1,o);
wind = -2:2
for i = 1:o
    si1 = i1(:,:,i);
    si1 = padthis(si1,3);
    stdv(i) = std(si1(:));
    meav(i) = mean(si1(:));
    [row,col] = find(si1 == maxv(i));
    sumv(i) = sum(sum(si1(row+wind, col + wind)));
    si1(row+wind, col + wind) = 0;
    resv(i) = sum(sum(si1));
end
thrsh = [meav;stdv;maxv;minv;resv;sumv];
% plot(meav+2*sqrt(meav),'r');
% hold on
% plot(meav,'o');
% plot(maxv,'y');
% plot(minv,'g');
% hold off
% figure
% plot(maxv./minv)
% figure
% plot(sumv)
% hold on
% plot(resv)
% hold off
% figure
% plot(sumv./resv)
