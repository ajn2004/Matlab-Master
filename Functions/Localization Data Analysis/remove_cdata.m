function cdata = remove_cdata(cdata,color, index)

field_arrays = fieldnames(cdata.(color)); % Get field array names
for i = 1:numel(field_arrays) % cycle through names
    if numel(cdata.(color).(field_arrays{i})) > 1 % avoid single entry names
        if numel( cdata.(color).(field_arrays{i})(1,:) ) == 6 % treat 2d arrays differently
            cdata.(color).(field_arrays{i})(index,:) = [];
        else
            cdata.(color).(field_arrays{i})(index) = []; % remove indices
        end
    end
end
