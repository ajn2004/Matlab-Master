function [ocoords] = astig_tilt(icords, cal)
% This function returns the corrected coordinates of an astigmatic system
% that has an inherent tilt as a function of Z. It requires the initial
% coordinates[x,y,z] (icords) the [x,y] tilts (tilts). it will return the
% corrected values xf and yf that are determined
zs = cal.tilt.zs;


xc = spline(zs(1:end-1) + mean(diff(zs))/2,cal.tilt.x,icords(:,3));
yc = spline(zs(1:end-1) + mean(diff(zs))/2,cal.tilt.y,icords(:,3));
xf = icords(:,1) - xc; 
yf = icords(:,2) - yc;
zf = icords(:,3)/cal.a(1);
ocoords = [xf,yf,zf];