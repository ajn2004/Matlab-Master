function ncoords = make_astigmatism_corrections(data,cal,q)
% A quick code to make the spline corrections to the astigmatism data.
% Data should be in 'pixel spatial units' with q = pixel/micron
% ncoords is returned in microns as the corrected data.
data(:,3) = data(:,3)/cal.a(1);
[xc, yc] = get_axial_tilt_spline(cal,data(:,3));
ncoords(:,1) = (data(:,1) - xc)*q;
ncoords(:,2) = (data(:,2) - yc)*q;
ncoords(:,3) = data(:,3)*q;
end