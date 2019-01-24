[m,n,o] = size(i1);
fms = 2;
for i = 1:round(o/3)
    im1(:,:,i) = mean(i1(:,:,1+(i-1)*fms:i*fms),3);
end