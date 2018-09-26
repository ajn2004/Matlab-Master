% Quick View
close all
i1 = readtiff()/33.33;
im2 = i1(:,:,2) - i1(:,:,1);
imagesc(im2.*(im2>0));