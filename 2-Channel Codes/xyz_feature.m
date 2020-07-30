function vector = xyz_feature(x,y,z)
x = x(:);
y=y(:);
z = z(:);
% vector = [x.^3, x.^2.*y, x.*y.^2, x.^2, y.^2, x.*y,x, y, x*0+1];
vector = [x.^2, y.^2, z.^2, x.*z, y.*z, x.*y, x, z, y, x*0+1];
vector = [x, z, y, x*0+1];
%vector = [x.^3, y.^3, x.^2.*y, x.*y.^2, x.^2, y.^2, x.*y,x, y, x*0+1];
    