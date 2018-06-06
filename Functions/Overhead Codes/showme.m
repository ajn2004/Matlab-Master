function showme(fname, fnum)
i1 = readtiff(fname,fnum,fnum+1);
imagesc(i1(:,:,1));
