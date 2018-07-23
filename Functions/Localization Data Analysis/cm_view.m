function cm_view(iloc)

[m,n,o] = size(iloc);
[X,Y] = meshgrid(1:m);
Xcm = sum(sum(X.*iloc));
Ycm = sum(sum(Y.*iloc));
xcm = Xcm/sum(iloc(:));
ycm = Ycm/sum(iloc(:));

imagesc(iloc);
hold on
plot(xcm,ycm,'k.')
hold off