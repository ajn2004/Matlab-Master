function i2 = gpu_conv(i1,kernel)
[m,n,o] = size(i1);
pixels = m*n*o; % number of pixels in image
i1 = double(i1); % ensure the image is of type double
maxpix = 224*228*1000; % empirical memory maximum before violating gpu
i2 = [];

if pixels > maxpix % if there are more pixels than the memory limit
    chunks = ceil(pixels/maxpix);  % break the process up into a number of chunks
    imchunks = round(o/chunks);
    for i = 1:chunks + 1
        im2 = [];
        try
           [im2] = gpu_conv_switch(i1(:,:,imchunks*(i-1)+1:imchunks*i),kernel);
        catch
            if imchunks*(i-1) +1 < o
                [im2] = gpu_conv_switch(i1(:,:,imchunks*(i-1)+1:end),kernel);
            end
        end
        i2 = cat(3,i2,im2);
    end
    
else
    [i2] = gpu_conv_switch(i1,kernel);
end
