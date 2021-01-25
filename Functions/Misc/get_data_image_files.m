function filtered_file_list = get_data_image_files(folder)
% A function to get a list of files that have data on them for
% localziation
cd(folder)
image_files = dir([folder,'*.tif']);
% Build cell of image files w/ data
    filtered_file_list = {};
    for l = 1:numel(image_files)
        if strfind(image_files(l).name, 'scan')
        else
            filtered_file_list{numel(filtered_file_list)+1,1} = image_files(l).name;
        end
    end
end