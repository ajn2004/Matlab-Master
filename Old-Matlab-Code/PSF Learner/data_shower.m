good = randi(numel(y)/2);

bad = randi(numel(y)/2)+numel(y)/2;
% create theta
for i = 1:7
im1(i,:) = theta(2+(i-1)*7:1+i*7);
img(i,:) =  X(good,2+(i-1)*7:1+i*7);
imb(i,:) = X(bad,2+(i-1)*7:1+i*7);
end

% subplot(1,3,1); imagesc(im1); colormap('gray'); title('Image of weight parameters'); axis image
subplot(1,2,1); imagesc(img); colormap('gray'); title('Image of passed tolerance molecule'); xlabel(['Confidence level is at ', num2str(100*sigmoid(X(good,:)*theta)), '%']); axis image
subplot(1,2,2); imagesc(imb); colormap('gray'); title('Image of failed tolerance molecule'); xlabel(['Confidence level is at ', num2str(100*sigmoid(X(bad,:)*theta)), '%']); axis image

figure
subplot(1,2,1); imagesc(img); colormap('gray');   axis image
subplot(1,2,2); imagesc(imb); colormap('gray');   axis image

