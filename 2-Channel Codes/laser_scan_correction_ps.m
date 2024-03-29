function t = laser_scan_correction_ps(file_list)

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
str_ind = strfind(ruler_name, 'd_');
str_end = strfind(ruler_name, '_r_');
step_size = str2num(ruler_name(str_ind+2:str_end-1))/1000;
if isempty(step_size)
    step_size = 0.02; % make default step size 20nm as that's what it set to default on the python gui 11/10/20 AJN
end


ruler_image = readtiff(ruler_name);
image = readtiff(image_name);
[m,n,o] = size(ruler_image);
channel_vertical_crop_stop = m;
sub_ruler_image = ruler_image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,:);
laser_images = image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,2:2:end);


indexes = get_index_from_ruler(sub_ruler_image,laser_images);
smoothed_indexes = conv(indexes,[0.2, 0.2, 0.2, 0.2, 0.2,],'same');
smoothed_indexes(1) = indexes(1);
smoothed_indexes(end) = indexes(end);
cal.z_drift = indexes; % we'll store the raw 
% Z correction for red
try
% z_correction_smoothed = spline((2:2:2*numel(indexes)), gausssmooth(indexes, 5,10),cdata.red.framenumber)*step_size;

z_correction_smoothed = spline((2:2:2*numel(indexes)), smoothed_indexes,cdata.red.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(indexes)), indexes,cdata.red.framenumber)*step_size;
cdata.red.zf_smoothed = cdata.red.zf - z_correction_smoothed;
cdata.red.zf_raw = cdata.red.zf - z_correction_raw;
catch
end
% Z correction for orange
try
% z_correction_smoothed = spline((2:2:2*numel(indexes)), smoothed_indexes,cdata.orange.framenumber)*step_size;
z_correction_smoothed = spline((2:2:2*numel(indexes)), smoothed_indexes,cdata.orange.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(indexes)), indexes,cdata.orange.framenumber)*step_size;
cdata.orange.zf_smoothed = cdata.orange.zf - z_correction_smoothed;
cdata.orange.zf_raw = cdata.orange.zf - z_correction_raw;
catch
end
% save([file_list{1}(1:end-4),'_sc.mat'],'cdata','cal','tol');
save([file_list{1}(1:end-4),'_old.mat'],'cdata','cal');
t = toc;
disp(['Corrected ', num2str(o),' frames in ', num2str(t), 's']);
catch lasterr
    waitforbuttonpress
    t=toc;
end
end

