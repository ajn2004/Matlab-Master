function red_indexes = get_index_from_curve(ruler_curve, laser_images, x_split)
% This function will take in a reference curve and laser bleedthrough
% images to establish the likely axial drift indices
[m,n,o] = size(laser_images);
summed_laser = sum(laser_images);
[m2] = size(ruler_curve);
clear laser_y laser_x
smoothed_difference_data = ruler_curve(20:180); % Grab the 'center' of the curve to avoid edge effects
[m2] = size(smoothed_difference_data);
% This will require review of every correction to ensure that the 'center'
% is in a response regime (i.e. no vertical section in curve v index)
for i = 1:o % Looping over every image to 
    laser_x = sum(summed_laser(1,1:x_split,i));
    laser_y = sum(summed_laser(1,x_split+1:end,i));
    laser_difference = (laser_x - laser_y)./(laser_x + laser_y);
    if abs(laser_difference - mean([min(smoothed_difference_data), max(smoothed_difference_data)])) < diff([min(smoothed_difference_data), max(smoothed_difference_data)])
        correction(i) = interp1(smoothed_difference_data,1:m2,laser_difference);
    else
        correction(i) = -1;
    end
end
indexes = find(correction == -1);
for i = 1:numel(indexes) % Set the 'correction' value in the event of a -1 (failed) event to be the previously found corrected position
    correction(indexes(i)) = correction(indexes(i)-1);
end
red_indexes = correction;
end