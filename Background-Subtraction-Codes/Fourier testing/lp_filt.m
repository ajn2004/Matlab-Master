function i2 = lp_filt(i1,coff)
[m,n,o] = size(i1);
[X,Y] = meshgrid((1:n) -n/2,(1:m)-m/2);
R = (X.^2 + Y.^2).^0.5;
G = exp(-R.^2/(2*coff^2));
G = G./(max(G(:)));
% [row, col] = find(R > coff);
for i = 1:o
    fi1 = fftshift(fft2(i1(:,:,i)));
    fi1 = fi1.*G;
    i2(:,:,i) = abs(ifft2(fftshift(fi1)));
end