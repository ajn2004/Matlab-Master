function [ruler_curve, x_split] = get_ruler_curve_from_images(sub_ruler_image)
% THis function will determine an 'ideal' separation on the ruler images to
% determine the comparator value (xl-xr)/(xl + xr)
[m, n, o] = size(sub_ruler_image);
summed_ruler_image = sum(sub_ruler_image);
mid_index = round(o/2);
% This loop will walk over columns summing all value before i and comparing
% them to the total sum overall columns
for i = 1:n
    photons(i) = sum(summed_ruler_image(1,1:i,mid_index))/sum(summed_ruler_image(1,:,mid_index));
end
photo_distance = abs(photons - 0.5); % max value of photons is 1, the value we want is the index where this curve is almost 0.5, here we're looking at the absolute distance to 0.5
% Grab the 'splitting' column position
x_split = find(photo_distance == min(photo_distance));
% Build the comparator curve using x_split as the divider of the columns
for i = 1:o
    xl(i,1) = sum(summed_ruler_image(1,1:x_split,i));
    xr(i,1) = sum(summed_ruler_image(1,x_split+1:end,i));
end
%Gaussian smooth the curve
x = -5:5;
gx = exp(-x.^2/(2*2.5^2));
gx = gx/sum(gx);
comparator_curve = (xl-xr)./(xl+xr);
% Pad the sides of the function to ensure good smoothing behavior
comparator = [ones(5,1)*comparator_curve(1);comparator_curve;ones(5,1)*comparator_curve(end)];
difference_data = conv(comparator,gx,'same');
%Grab only original curve indices
ruler_curve = difference_data(6:end-5);
end
