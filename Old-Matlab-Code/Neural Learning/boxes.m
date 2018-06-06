function boxes(cents,width,color)
% This function will draw boxes onto an open image centered at cents with
% smallest radius = width and with color = color This assumes an image is
% open
hold on
for i = 1:numel(cents(:,1));
    x0 = cents(i,1) - width;
    x1 = cents(i,1) + width;
    y0 = cents(i,2) - width;
    y1 = cents(i,2) + width;
    plot([x0, x0, x1, x1, x0],[y0, y1,y1,y0, y0],color);
end
hold off
end