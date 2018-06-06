function circles(x,y,r,c)
% this function will take in 3 vectors of x, y, and r of equal size and
% plot circles of radius r(i) with centers (x(i),y(i))
% c will determine the color of the circle used

% AJN 9/22/15
hold on
th = 0:pi/100:2*pi;
for i = 1:numel(x)
    xunit = r(i) * cos(th) + x(i);
    yunit = r(i) * sin(th) + y(i);
    plot(xunit, yunit,c, 'MarkerSize', 20);
   
end
hold off
end