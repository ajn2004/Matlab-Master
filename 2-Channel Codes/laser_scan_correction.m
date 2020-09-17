function t = laser_scan_correction(file_list)

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
data_name = [image_name(1:end-4),'_dast_tol.mat'];
ruler_image = readtiff(ruler_name);
image = readtiff(image_name);
[m,n,o] = size(ruler_image)
channel_vertical_crop_stop = m;
sub_ruler_image = ruler_image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,:);
laser_images = image(channel_vertical_crop_start:channel_vertical_crop_stop,channel_split:channel_split_end,2:2:end);


indexes = get_index_from_ruler(sub_ruler_image,laser_images);

% Z correction for red
z_correction_smoothed = spline((2:2:2*numel(indexes)), gausssmooth(indexes, 5,10),cdata.red.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(indexes)), indexes,cdata.red.framenumber)*step_size;
cdata.red.zf_smoothed = cdata.red.zf - z_correction_smoothed;
cdata.red.zf_raw = cdata.red.zf - z_correction_raw;

% Z correction for orange
z_correction_smoothed = spline((2:2:2*numel(indexes)), gausssmooth(indexes, 5,10),cdata.orange.framenumber)*step_size;
z_correction_raw = spline((2:2:2*numel(indexes)), indexes,cdata.orange.framenumber)*step_size;
cdata.orange.zf_smoothed = cdata.orange.zf - z_correction_smoothed;
cdata.orange.zf_raw = cdata.orange.zf - z_correction_raw;

save(['scan_corrected\',file_list{1}(1:end-4),'_sc.mat'],'cdata','cal','tol');
t = toc;
disp(['Corrected ', num2str(o),' frames in ', num2str(t), 's']);
catch lasterr
    waitforbuttonpress
    t=toc;
end
end

function indexes = get_index_from_ruler(ruler, image_stack)
    [m,n,o] = size(image_stack);
    [m,n,r] = size(ruler);
    indexes = zeros(o,1);
    for i = 1:o % loop over images
        if round(100*i/o) == 100*i/o
            clc;
            disp(['Indexed ', num2str(100*i/o), '% of file']);
        end
        C = [];
        for j = 1:r
           C(j) = corr2(ruler(:,:,j), image_stack(:,:,i));  
        end
        
        indexes(i) = find(C == max(C));

    end
end