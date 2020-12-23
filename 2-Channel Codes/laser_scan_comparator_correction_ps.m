function t = laser_scan_correction_ps(file_list)
% This function is to attempt to analyze the laser bleedthrough image to
% determine the necessary axial drift correction for a data set
error = 0; % error will be used as a readout in the instance somethign goes wrong
try
    tic
% These values should be updated based on calibrations and changes to
% detection path
channel_split = 200; % X position split in pixels
channel_split_end = 375; % end of channel cropping
channel_vertical_crop_start = 15;



% If I code the rest of this correctly the user doesn't need to interact
% with the code past this point.

ruler_name = file_list{3};
image_name = file_list{2};
load(file_list{1});
error = 1; % Error 0 indicates a problem w/ loading the localized data set
str_ind = strfind(ruler_name, 'd_');
str_end = strfind(ruler_name, '_r_');
step_size = str2num(ruler_name(str_ind+2:str_end-1))/1000;
if isempty(step_size)
    step_size = 0.02; % make default step size 20nm as that's what it set to default on the python gui 11/10/20 AJN
end

% error 1 indicates issue w/ reading the ruler name
ruler_image = readtiff(ruler_name);
error = 2; % Error 2 indicates issue w/ loading the raw images
image = readtiff(image_name);
error = 3; % Error 3 indicates issue w/ getting ruler from curve
[m,n,o] = size(ruler_image);
channel_vertical_crop_stop = m;
% Simple cropping step
sub_ruler_image = ruler_image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,:);
laser_images = image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,2:2:end);
[ruler_curve, x_split] = get_ruler_curve_from_images(sub_ruler_image); % Build the (x-y)/(x+y) reference curve
error = 4; % Indicates issue w/ getting raw indicies from ruler curve
% indexes = get_index_from_ruler(sub_ruler_image,laser_images);
red_indexes = get_index_from_curve(ruler_curve, laser_images,x_split); % Build the indexed correction
smooth_red = [ones(1,5)*red_indexes(1),red_indexes,ones(1,5)*red_indexes(end)];
smoothed_indexes = conv(smooth_red,[0.2, 0.2, 0.2, 0.2, 0.2,],'same');
smoothed_red_indexes = smoothed_indexes(6:end-5);
cal.z_drift = red_indexes; % we'll store the raw 
% Z correction for red
try
% z_correction_smoothed = spline((2:2:2*numel(indexes)), gausssmooth(indexes, 5,10),cdata.red.framenumber)*step_size;

z_correction_smoothed = spline((2:2:2*numel(red_indexes)), smoothed_red_indexes,cdata.red.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(red_indexes)), red_indexes,cdata.red.framenumber)*step_size;
cdata.red.zf_smoothed = cdata.red.zf - z_correction_smoothed;
cdata.red.zf_raw = cdata.red.zf - z_correction_raw;
catch
end
% Z correction for orange
try
% z_correction_smoothed = spline((2:2:2*numel(indexes)), smoothed_indexes,cdata.orange.framenumber)*step_size;
z_correction_smoothed = spline((2:2:2*numel(red_indexes)), smoothed_red_indexes,cdata.orange.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(red_indexes)), red_indexes,cdata.orange.framenumber)*step_size;
cdata.orange.zf_smoothed = cdata.orange.zf - z_correction_smoothed;
cdata.orange.zf_raw = cdata.orange.zf - z_correction_raw;
catch
end
% save([file_list{1}(1:end-4),'_sc.mat'],'cdata','cal','tol');
save([file_list{1}(1:end-4),'_sc.mat'],'cdata','cal');
t = toc;
disp(['Corrected ', num2str(o),' frames in ', num2str(t), 's']);
catch lasterr
    disp('Issue in laser scan correction')
    switch error
        case 0
            disp(['Problem likely exists in data set named ', file_list{1}])
        case 1
            disp(['Problem likely exists in image named ', ruler_name])
        case 2
            disp(['Problem likely exists in image named ', image_name])
        case 3
            disp('Problem is likely an issue in building the reference curve')
        case 4
            disp('Problem likely is an issue w/ getting indices from reference curve')
        otherwise
            disp('Time for some thrilling debugging')
    end
    waitforbuttonpress
    t=toc;
end
end

