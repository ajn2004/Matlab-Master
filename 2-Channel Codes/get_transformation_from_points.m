function [o2rx, o2ry] = get_transformation_from_points(data)
% expecting data to come as [xo, yo,xr,yr]
xo = data(:,1);
yo = data(:,2);
xr = data(:,3);
yr = data(:,4);
MX = xy_feature(xo,yo);
MY = xy_feature(xo,yo);
MX = [MX,xr];
MY = [MY,yr];
MX = rref(MX);
MY = rref(MY);
o2rx = MX(:,end);
o2ry = MY(:,end);