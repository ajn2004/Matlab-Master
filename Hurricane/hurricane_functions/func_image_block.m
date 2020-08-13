function i2 = func_image_block(i1,split,toggle)
[m,n,o] = size(i1);
i2 = i1*0;
oblock = ones(m,n);
rblock = oblock;

oblock(:,split:end) = 0;
rblock(:,1:split) = 0;

if toggle == 1
    for i = 1:2:o
        i2(:,:,i) = i1(:,:,i).*rblock;
    end
    
    for i = 2:2:o
        i2(:,:,i) = i1(:,:,i).*oblock;
    end
    
else
    for i = 2:2:o
        i2(:,:,i) = i1(:,:,i).*rblock;
    end
    
    for i = 1:2:o
        i2(:,:,i) = i1(:,:,i).*oblock;
    end
    
end