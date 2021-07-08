function ncoords = make_astigmatism_corrections(data,cal,q)
% We perform the corrections in 'pixels' before converting to microns
% The calibration software dealt directly with fitted localizations, so we
% carried that into this software. 
data(:,3) = (data(:,3) - cal.a(2))/cal.a(1); % Subtracting off cal.a(2) rectifies 3D offset between channels. 
% The calibration software deals with the magnification corrected
% coordinates
[xc, yc] = get_axial_tilt_spline(cal,data(:,3)); % get spline corrections

% Make corrections and output result in microns
ncoords(:,1) = (data(:,1) - xc)*q;
ncoords(:,2) = (data(:,2) - yc)*q;
ncoords(:,3) = data(:,3)*q;
% end