function set_scale(i1,q,x)
% sets the scale to reflect physical size
% i1 is the image, q is the pixel to um variable, x is the number of ticks
% to have

[m,n,o] = size(i1);
xtix = floor(n/x);
ytix = floor(m/x);
clims = [min(i1(:)),13000];

for i = 1:o
%     figure
    imagesc(i1(:,:,i),clims);
    set(gca,'Xtick',0:xtix:n);
    set(gca,'XtickLabel',(0:xtix:n)*q);
    set(gca,'Ytick',0:ytix:m);
    set(gca,'YtickLabel',(0:ytix:m)*q);
    title(['Image number ',num2str(i)]);
    xlabel('Position(um)');
    ylabel('Position(um)');
    colormap('jet');
end    
axis image