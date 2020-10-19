function i2 = func_image_block(i1,split,toggle)
% Quickly block out regions of data that we don't care about.
[m,n,o] = size(i1);
i2 = i1*0;
oblock = ones(m,n);
rblock = oblock;

oblock(:,split:end) = 0;
rblock(:,1:split) = 0;

if toggle == 1 % 1 for having red on the even frames and orange on the odd
    for i = 1:2:o
        i2(:,:,i) = i1(:,:,i).*rblock;
    end
    
    for i = 2:2:o
        i2(:,:,i) = i1(:,:,i).*oblock;
    end
    
else % 2 for having orange on the even and red on the odd frames
    for i = 2:2:o
        i2(:,:,i) = i1(:,:,i).*rblock;
    end
    
    for i = 1:2:o
        i2(:,:,i) = i1(:,:,i).*oblock;
    end
    
end