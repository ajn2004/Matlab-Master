function [corrected_data] = channel_correct(data)

% corrected_data = data*0; % initialize data
try
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat','split','o2rx','o2ry');
catch
    load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat','split','o2rx','o2ry');
end
xf = data(:,1);
yf = data(:,2);

% vector = [xfr.^2, yfr.^2, xfr, yfr,xfr.*yfr, xfr*0+1];
vector = xy_feature(xf,yf);
xfo = o2rx.'*vector.';
yfo = o2ry.'*vector.';

corrected_data = [xfo(:),yfo(:)];
end
