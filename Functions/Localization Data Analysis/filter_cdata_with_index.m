function cdata1 = filter_cdata_with_index(cdata,index)
cdata1 = cdata;

colors = {'red','orange'};
for i  = 1:2 % loop over colors
    fieldnames = fields(cdata.(colors{i})); % Get relevant field names for the object
    for j  = 1:2 % First 2 fields are 2D arrays
        cdata1.(colors{i}).(fieldnames{j}) = cdata.(colors{i}).(fieldnames{j})(index.(colors{i}),:);
    end
    for j = 3:numel(fieldnames) % remaining fields are all 1D
        cdata1.(colors{i}).(fieldnames{j}) = cdata.(colors{i}).(fieldnames{j})(index.(colors{i}));
    end
end