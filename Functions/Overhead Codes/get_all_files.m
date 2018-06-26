function [files] = get_all_files(ext,fpath)

ll = genpath(fpath); % generate a list of all subfolders
ind = strfind(ll,';'); % list is separated by semicolons
subs{1} = ll(1:ind(1)-1); % subs will be the variable of all subfolders
for i = 2:numel(ind)
    subs{i} = ll(ind(i-1)+1:ind(i)-1);
end

files.name = [];
files.folder = [];
files.date = [];
files.bytes = [];
files.isdir = [];
files.datenum = [];
for i = 1:numel(subs) %build files field variable by looping through sub folders
    flist =[];
    flist = dir([subs{i},'\*',ext]); % look for files with desired extension
    if ~isempty(flist)
        for j = 1:numel(flist)
            ind = numel(files) +1;
            files(ind) = flist(j);
        end
    end
end
files(1) = [];