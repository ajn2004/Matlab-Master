function vector = xy_feature(x,y)
x = x(:);
y=y(:);
% vector = [x.^3, x.^2.*y, x.*y.^2, x.^2, y.^2, x.*y,x, y, x*0+1];
<<<<<<< HEAD
% vector = [x.^2, y.^2, x.*y,x, y, x*0+1];
=======
vector = [x.^2, y.^2, x.*y,x, y, x*0+1];
>>>>>>> master
vector = [x.^3, y.^3, x.^2.*y, x.*y.^2, x.^2, y.^2, x.*y,x, y, x*0+1];
    