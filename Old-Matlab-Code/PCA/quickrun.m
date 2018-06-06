i1 = readtiff('storm.tif',200,203);
[m,n,o] = size(i1);
y = filtcut(i1,3,o);