function [xc,yc] = get_axial_tilt_spline(z_cal,z0)
% Looking through total calibrations shows that these parameters are
% expected in 'pixel space'
xc = spline(z_cal.z_tilt,z_cal.x_tilt,z0);
yc = spline(z_cal.z_tilt,z_cal.y_tilt,z0);
xc = xc(:);
yc = yc(:);