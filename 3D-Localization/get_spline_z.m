function zf_um = get_spline_z(sigma_x, sigma_y, z_cal)
% this function returns Z in microns from sigma widths in pixels
    x_spread = min(z_cal.z0s):0.001:max(z_cal.z0s); % Curve span in microns
    x_spread = -2:0.001:2;
    x_curve = spline(z_cal.z0s,gausssmooth(z_cal.sx,2,5),x_spread);
    y_curve = spline(z_cal.z0s,gausssmooth(z_cal.sy,2,5),x_spread); 
    zf_um = sigma_x*0;
    for i = 1:numel(sigma_x)
%     D = ((sigma_x(i).^0.5-x_curve.^0.5).^2 + (sigma_y(i).^0.5-y_curve.^0.5).^2).^0.5;
%     D = ((sigma_x(i)^2-x_curve.^2).^2 + (sigma_y(i)^2-y_curve.^2).^2).^0.5;
    D = (sigma_y(i)/sigma_x(i)-y_curve./x_curve).^2;
%     D = abs(x_curve./y_curve - sigma_x(i)/sigma_y(i));
    ind = find(D == min(D), 1);
    min_d(i) = D(ind);
    if D(ind) < 1 && (ind ~= 1 || ind ~= numel(x_spread))
%         plot(x_spread,D)
%         hold on
        try
            zf_um(i) = x_spread(ind(1));  % emperically found to be the axial zoom factor by fitting a line between board movements and position measurements
        catch
            zf_um(i) = -5;
        end
    else 
        zf_um(i) = -5;
    end
    
    end
disp(['Of the values entered, ' num2str(100*sum(zf_um==-5)/numel(zf_um)),'% were rejected']);