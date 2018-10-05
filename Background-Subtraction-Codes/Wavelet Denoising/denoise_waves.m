function i2 = denoise_waves(i1,k)
% This function will denoise all images in i1 to level k and return an
% array i2 of same size i1 w/ denoised values
[m,n,o]= size(i1);
for i = 1:o
[W, I] = get_waves(i1(:,:,i),k);
W1 =thresh_waves(W);

i2(:,:,i) = sum(W1,3)+I(:,:,end);
end