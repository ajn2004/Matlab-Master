function cdata = clean_cdata(cdata)
% A quick script to remove localizations whose fitting parameters == 0
    colors = {'red','orange'};
    
    for j = 1:numel(colors)
        try
        remove_index = (cdata.(colors{j}).xf == 0 & cdata.(colors{j}).yf == 0) | cdata.(colors{j}).zf < 0;
        disp(sum(remove_index))
        field_arrays = fieldnames(cdata.(colors{j})); % Get field array names
        for i = 1:numel(field_arrays) % cycle through names
            if numel(cdata.(colors{j}).(field_arrays{i})) > 1 % avoid single entry names
                if numel( cdata.(colors{j}).(field_arrays{i})(1,:) ) == 6 % treat 2d arrays differently
                    cdata.(colors{j}).(field_arrays{i})(remove_index,:) = [];
                else
                    cdata.(colors{j}).(field_arrays{i})(remove_index) = []; % remove indices
                end
            end
        end
        catch
        end
    end
end