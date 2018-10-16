function [ind, mc] = cross3d(i1,i2)
% A function to cross correlate 3D scans to get a more robust measurement
% of the defocusing curves of my system
% This function convolve i2 w/ i1 and return the index of the highest
% 3-index correlation AJN 9/19/2018

ind = [];
[m,n,~] = size(i1);

for i = 1:m % cross correlate over all rows and cols
    for j = 1:n
        i2s = i2(i,j,:);
        i33s = i1(i,j,:);
        C(i,j,:) = conv(i33s(:),i2s(:));
    end
end

mc = mean(mean(C));
ind = find(mc == max(mc)); % return index of largest <correlation>