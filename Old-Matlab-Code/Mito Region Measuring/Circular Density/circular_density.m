function [rs, ds] = circular_density(xc,yc, radius, i1, fact)
% This script will measure the density of a region as a function of radius
% for an image given It will measure from 0 to fact*radius and return the
% radii looked at and the value of the density at that point Density is
% measured by the sum of pixel values divided by the area covered by those
% pixels
%
% AJN 2/19/16

%preallocate
rs = 0:0.5:100;
ds = rs*0;

for i = 1:numel(rs)
    r = rs(i);
    xs = round(xc-r):round(xc+r);
    ys = round(yc-r):round(yc+r);
    try
    reg = i1(ys,xs);
    mass = sum(reg(:));
    [m,n] = size(reg);
    areas = pi*r*r;
    ds(i) = mass/areas;
    catch lsterr
    end
%     waitforbuttonpress
end