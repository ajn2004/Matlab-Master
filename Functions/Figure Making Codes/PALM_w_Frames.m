%% PALM W FRAMES
% clearvars; close all; clc
close all
% i1 = readtiff();
[m,n,o] = size(im1);
% im1 = roball(im1,6,6);
pixw = 5;
for i = 1:o
    ind = find(framenum_all == i);
    subplot(1,2,1);
    imagesc(im1(:,:,i));
    hold on
    draw_boxes([fits(ind,1),fits(ind,2)],pixw);
    axis image
    hold off
    subplot(1,2,2);
    if ~isempty(ind)
    plot(xf_fixed(1:ind(end))*q,-q*yf_fixed(1:ind(end)),'.b');
    hold on
    plot(q*xf_fixed(ind),-q*yf_fixed(ind),'.r');
    axis equal
    
    xlabel('Position in um');
    ylabel('Position in um');
    end
    drawnow
    M(i) = getframe(gcf);
end