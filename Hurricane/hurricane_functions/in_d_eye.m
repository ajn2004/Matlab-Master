function in_d_eye(iprod, dps, pixw)
% A quick function to draw boxes around rendered data as an example
[~, ~, o] = size(iprod);
for i = 1:o
%     subplot(1,2,1)
%     imagesc(i1(:,:,i))
%     axis image
%     title('Raw Data');
%     subplot(1,2,2);
    imagesc(iprod(:,:,i))
    [row,col] =find(dps(:,:,i) == 1);
    draw_boxes([col,row],pixw);
    colormap('gray')
    axis image
    title('Background Subtracted')
    drawnow
    M(i) = getframe(gcf);
end

movie2gif(M,'frames.gif','DelayTime',0.03,'LoopCount',Inf);
