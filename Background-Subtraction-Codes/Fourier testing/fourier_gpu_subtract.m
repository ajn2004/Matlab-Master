function iprod = fourier_gpu_subtract(i1,sp)
[m,n,o] = size(i1);
i2 = permute(i1,[3,2,1]);
clear i1
i3 = reshape(i2,o,m*n);
clear i2
fi = abs(cgpufourier2(i3,sp))/o;
m3 = mean(i3(:));
mi = mean(fi(:));
clear i3
fi = reshape(fi,o,n,m);
iprod = permute(fi,[3,2,1]);
clear fi rat mi m3
end
