function d = dist_meas();
[x,y] = ginput(2);
d = ((x(1)-x(2))^2 + (y(1)-y(2))^2)^0.5;
disp(d);