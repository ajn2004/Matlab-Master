function file_list_out = rename_problem_files(file_list_in)
% This function will look at the name of a data file and rename it so that
% there should be no issue w/ subsquent 'perfect storm' analysis
file_list_out = file_list_in;
for i = 1:numel(file_list_in)
    % don't rename *scan* files, those are usually named responsibly
    try
        % we want the file name to end in _data.tif
        ind1 = strfind(files(i).name,'_r_');
        if ~isempty(ind1) % If an index was found it's because we've incorrectly labeled a data file
            if isempty(strfind(files(i).name,'scan'))
                new_name = [files(i).name(1:ind1),'data.tif'];
                movefile(files(i).name, new_name)
                file_list_out(i).name = new_name;
            else
                new_name = [files(i).name(1:ind1),'scan.tif'];
                movefile(files(i).name, new_name)
                file_list_out(i).name = new_name;
            end
        end
    catch
    end
    
end
end