function im2 = gpu_waves(im1)
% This code is to process im1 in a wavelet denoising using gpu acceleration
% The waveframe we're after is frame 2 so we need to blur 2x then take the
% difference of the two blurrs and return that

im2 = im1*0; % preallocate
baselet = [1/16, 1/4, 3/8, 1/4, 1/16];

i = 1;
wvlt = [baselet(1), zeros(1,2^(i-1)-1),baselet(2), zeros(1,2^(i-1)-1),baselet(3), zeros(1,2^(i-1)-1),baselet(4), zeros(1,2^(i-1)-1),baselet(5)]; % construct wavelet
wave = wvlt.*wvlt.';
i1 = gpu_conv(im1,wave);
disp('Wavelet: Conv step 1 done')
i = 2;
wvlt = [baselet(1), zeros(1,2^(i-1)-1),baselet(2), zeros(1,2^(i-1)-1),baselet(3), zeros(1,2^(i-1)-1),baselet(4), zeros(1,2^(i-1)-1),baselet(5)]; % construct wavelet
wave = wvlt.*wvlt.';
i2 = gpu_conv(i1,wave);
disp('Wavelet: Conv step 2 done')
im2 = split_sub(i1,i2);