function folders = get_all_analysis_folders(current_folder)
% Let's only grab analysis folders since those will necessarily come
% from data folders
f = genpath(current_folder);
folders = {};
index = strfind(f,';');
for i = 1:numel(index)-1
    an_index = strfind(f(index(i):index(i+1)),'Analysis;');
    if ~isempty(an_index)
        folders{numel(folders)+1,1} = [f(index(i)+1:index(i+1)-1),'\'];
    end
end
end