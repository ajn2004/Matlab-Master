function set_scale_ax(i1,q,x,ax)
% sets the scale to reflect physical size
% i1 is the image, q is the pixel to um variable, x is the number of ticks
% to have

    
[m,n,o] = size(i1);
xtix = floor(n/x);
ytix = floor(m/x);


for i = 1:o

    imagesc(ax,i1(:,:,i));
    set(gca,'Xtick',0:xtix:n);
    set(gca,'XtickLabel',(0:xtix:n)*q);
    set(gca,'Ytick',0:ytix:m);
    set(gca,'YtickLabel',(0:ytix:m)*q);
    title(['Image number ',num2str(i)]);
    xlabel('Position(um)');
    ylabel('Position(um)');
    colormap('jet');
end    
    