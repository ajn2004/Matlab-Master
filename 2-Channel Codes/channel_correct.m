function [corrected_data] = channel_correct(data)

% corrected_data = data*0; % initialize data
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat','split','o2rx','o2ry');
xf = data(:,1);
yf = data(:,2);
ind = xf<split;
xfr = xf(ind);
yfr = yf(ind);
ind = logical(1-ind);
xfo = xf(ind);
yfo = yf(ind);

% vector = [xfr.^2, yfr.^2, xfr, yfr,xfr.*yfr, xfr*0+1];
vector = xy_feature(xfo,yfo);
xfo = o2rx.'*vector.';
yfo = o2ry.'*vector.';


cdata.xf{1} = xfr;
cdata.xf{2} = xfo;
cdata.yf{1} = yfr;
cdata.yf{2} = yfo;
corrected_data = cdata;
