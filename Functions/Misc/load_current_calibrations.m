function [cal] = load_current_calibrations()
%load_current_zcalibration Summary of this function goes here
%   Detailed explanation goes here
try
    load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
    load('C:\Users\andre\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat')
catch
    
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat')
end
cal.o2rx = o2rx;
cal.o2ry = o2ry;
cal.split = split;
end

