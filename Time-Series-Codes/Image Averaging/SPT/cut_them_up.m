function [i2s, framenum, beta0] = cut_them_up(i1, points, psize)
% Initialize variables
i2s =[];
framenum =[];
beta0 = [];
count = 1;
[xgrid, ygrid] = meshgrid(-psize:psize,-psize:psize);

% loop over all frames for the region
for j = 1:numel(i1(1,1,:))

        i2s(:,:,count) = i1(points(1,2) - psize : points(1,2) + psize, ...
            points(1,1) - psize : points(1,1) + psize,j);
        framenum(count) = j;
        xcm = sum(sum(xgrid.*i2s(:,:,count)))/numel(xgrid(:));
        ycm = sum(sum(ygrid.*i2s(:,:,count)))/numel(xgrid(:));
        beta0(count,1) = xcm;
        beta0(count,2) = ycm;
        count = count+1;
    end



end