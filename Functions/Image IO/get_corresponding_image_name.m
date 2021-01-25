function corresponding_image_name = get_corresponding_image_name(filesname,data_image_files)
    for l = 1:numel(data_image_files)
        if strfind(data_image_files{l},filesname(1:8))
            corresponding_image_name = filesname;
        end
    end
    end