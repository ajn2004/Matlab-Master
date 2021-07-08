function combine_cdata_folder(folder)

files = dir([folder,'*dast.mat']);
% Build file sequences
% our naming convention will never contain X2 in the base name of the file
% We cna use regular expressions to extract this information

for i = 1:numel(files)
index = regexp(files(i).name,'X[1-9]');
    if ~isempty(index) % if index is empty, it's the first file of a set
        basename = [files(i).name(1:index-2),'_dast.mat'];
        for j = i:numel(files)
            if strcmp(files(j).name,basename)
                load([folder,files(j).name]);
                cdata_one = cdata;
                load([folder,files(i).name]);
                cdata_two = cdata;
                clearvars cdata;
                cdata = combine_cdata(cdata_one,cdata_two);
                save([folder,files(j).name(1:end-4),'_combo.mat']);
                break
            end
        end
    end


end
end

