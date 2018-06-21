function [snr, sig, nois] = sig2noi(i1,i2,aves)

[m,n,o] = size(i1);

for i = 1:o
    sig(i) = sum(sum(i1(:,:,i)));
    bkgn = i2(:,:,(i-1)*aves +1 : i*aves);
    nois(i) = std(bkgn(:));
end
snr = mean(sig./nois);