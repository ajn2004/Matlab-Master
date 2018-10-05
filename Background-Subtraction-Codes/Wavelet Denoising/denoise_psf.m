function i2 = denoise_psf(i1,ts)
% This function will denoise all images in i1 to level k and return an
% array i2 of same size i1 w/ denoised values
[m,n,o]= size(i1);
for i = 1:o
[W, I] = get_waves(i1(:,:,i),3);
% Thresholding
[m,n,o] = size(W);
[fits] = fit_hist_gauss(W(:,:,1));
lvls = o;
% Thresholding the components
% Thresholding is done to the W matrix
for k = 1:lvls
    ii3 = W(:,:,k);
    thrsh = ts*abs(fits(3));
    W1(:,:,k) = W(:,:,k).*(W(:,:,k)>thrsh);
end
i2(:,:,i) = W1(:,:,2);
end