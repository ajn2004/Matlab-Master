function draw_circ(x,y,r)
x1 = linspace(-r,r,1000);
y1 = sqrt(r^2 - (x1).^2) - y;
hold on
plot(x1+x,y1, 'r');
plot(x1+x,-y1,'r');
hold off