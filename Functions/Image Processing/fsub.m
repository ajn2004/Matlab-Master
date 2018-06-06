function i3 = fsub(i1, sp)
%FSUB Using the fourier method to subtract background.
%   I2 = FSUB(IMAGE, SETPOINT) takes in a image series specified by IMAGE.
%   This algorithm assumes the 3rd index as the temporal index. The
%   function will the apply a high pass filter at index SETPOINT which can
%   be converted to frequency by dividing by the spacing in your data set.

% Find size of image
[m, n, p] = size(i1);
% i2 = zeros(m,n,p);
% Preallocation and fft
nch = 10;
chs = floor(m/nch);
lef = mod(m,nch);
for i = 1:nch
i2(1:chs,:,:) = fft(i1((i-1)*chs +1:i*chs,:,:),p,3);
% clear i1
%% High Pass Filter
 i2(:,:,1:sp) = i2(:,:,1:sp)*0;
 i2(:,:,end-sp+1:end) = i2(:,:,end-sp+1:end)*0;
 % Reconstruction of image space
 i3((i-1)*chs +1:i*chs,:,:) = abs(ifft(i2,p,3));
 clear i2
end
if lef > 0
    i2(1:lef,:,:) = fft(i1(i*chs+1:end,:,:),p,3);
%% High Pass Filter
 i2(:,:,1:sp) = i2(:,:,1:sp)*0;
 i2(:,:,end-sp+1:end) = i2(:,:,end-sp+1:end)*0;
 % Reconstruction of image space
 i3(i*chs+1:m,:,:) = abs(ifft(i2,p,3));
 clear i2
end
