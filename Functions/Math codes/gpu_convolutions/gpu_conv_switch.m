function i2 = gpu_conv(i1,k1)
% This function allows for the gpu convolution of image i1 with kernel k1

[m,n,o] = size(i1);
[mm, nn] = size(k1);

if mm == nn
    switch mm
        case 3
            i2 = gpu_conv_3(i1,k1);
        case 5
            i2 = gpu_conv_5(i1,k1);
        case 7
            i2 = gpu_conv_7(i1,k1);
        case 9
            i2 = gpu_conv_9(i1,k1);
        case 11
            i2 = gpu_conv_11(i1,k1);
        case 13
            i2 = gpu_conv_13(i1,k1);
        case 15
            i2 = gpu_conv_15(i1,k1);            
        case 17
            i2 = gpu_conv_17(i1,k1);
        case 19
            i2 = gpu_conv_19(i1,k1);
        case 21
            i2 = gpu_conv_21(i1,k1);
        case 23
            i2 = gpu_conv_23(i1,k1);
    end
end