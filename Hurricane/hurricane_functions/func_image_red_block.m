function i2 = func_image_red_block(i1,split)
[m,n,o] = size(i1);
i2 = i1*0;
oblock = ones(m,n);
rblock = oblock;

oblock(:,split:end) = 0;
rblock(:,1:split) = 0;

for i = 1:o
    i2(:,:,i) = i1(:,:,i).*rblock;
end
