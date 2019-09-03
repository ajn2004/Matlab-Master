function i1 = func_3D_hist(coord,bin, radius)
for i = 1:3
    coord(:,i) = coord(:,i) - min(coord(:,i));
end
xf = coord(:,1);
yf = coord(:,2);
zf = coord(:,3);

sbin = bin/1000;
xl = [min(xf), max(xf)];
yl = [min(yf), max(yf)];
zl = [min(zf), max(zf)];
m = (floor(diff(yl)/sbin)+1);
n = (floor(diff(xl)/sbin)+1);
o = (floor(diff(zl)/sbin)+1);
% i1 = int32(zeros(m,n,o)); % Create 3D histogram
xf = (xf - xl(1))*n/xl(2);
yf = (yf - yl(1))*m/yl(2);
zf = (zf - zl(1))*o/zl(2);
rads = radius*n/xl(2);
% Populate the 3D Histogram
i1 = zeros(m,n,o);

for i = 1:m
    for j = 1:n
        for k = 1:o
            ind = (xf-j).^2 + (yf - i).^2 + (zf - k.^2) <= rads^2;
            i1(i,j,k) = sum(ind);
        end
    end
end

 

