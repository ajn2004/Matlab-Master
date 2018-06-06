function h = circles(x,y,r,c)
% this function will take in 3 vectors of x, y, and r of equal size and
% plot circles of radius r(i) with centers (x(i),y(i))
hold on
th = 0:pi/100:2*pi;
for i = 1:numel(x)
    xunit = r(i) * cos(th) + x(i);
    yunit = r(i) * sin(th) + y(i);
    h = plot(xunit, yunit,c);
    hold off
end
end