function zf_um = get_spline_z(sigma_x, sigma_y, z_cal)
% this function returns Z in microns from sigma widths in pixels
    x_spread = -0.5:0.001:0.5; % Curve span in microns
    x_curve = spline(z_cal.z0s,z_cal.sx,x_spread);
    y_curve = spline(z_cal.z0s,z_cal.sy,x_spread); 
    zf_um = sigma_x*0;
    for i = 1:numel(sigma_x)
    D = ((sigma_x(i).^0.5-x_curve.^0.5).^2 + (sigma_y(i).^0.5-y_curve.^0.5).^2).^0.5;
%     D = ((sigx(i)-sx).^2 + (sigy(i)-sy).^2).^0.5;
%     D = ((sigx(i).^0.5-sx.^0.5).^2 + (sigy(i).^0.5-sy.^0.5).^2);
    ind = find(D == min(D), 1);
    try
    zf_um(i,1) = x_spread(ind(1));  % emperically found to be the axial zoom factor by fitting a line between board movements and position measurements
    catch
        zf_um(i,1) = -300000;
    end
    end
end