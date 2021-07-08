function cdata_out = combine_cdata(cdata_one,cdata_two)
%a function to combine multiple cdata structures (while removing PSF)
fnames = fieldnames(cdata_one.red);
if strcmp(fnames{5} , 'psfs') % check to see if PSF field name exists
    fnames(5) = [];% remove it if found
end
% Prepare output data structure
cdata_out.red = [];
cdata_out.orange = [];
% loop through colors
for j = {'orange','red'}
    % Loop through field names
    for i = 1:numel(fnames)
        if i == 4 % Every field can be concatenated except framenumber, which needs to have an offset for the second structure (structure two starts after structure one)
            cdata_out.(j{1}).(fnames{i}) = [cdata_one.(j{1}).(fnames{i})(:);max(cdata_one.(j{1}).(fnames{i})) + cdata_two.(j{1}).(fnames{i})(:)];
        else
            cdata_out.(j{1}).(fnames{i}) = [cdata_one.(j{1}).(fnames{i}); cdata_two.(j{1}).(fnames{i})];
        end
    end
end
end

