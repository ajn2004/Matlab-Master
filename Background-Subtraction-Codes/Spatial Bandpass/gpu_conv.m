function i2 = gpu_conv(i1,k1)
% This function allows for the gpu convolution of image i1 with kernel k1

[m,n,o] = size(i1);
[mm, nn, oo] = size(k1);

if mm == nn && mm < 15
    i2 = image_conv(single(i1),single(k1));
end