clearvars
try
base_folder = 'G:\Dropbox\Andrew & Tim shared folder\Analyzed Data\';

data_folder = pwd;
index = strfind(data_folder, 'Data');
endex = strfind(data_folder, 'Analysis');

new_save_folder = [base_folder, data_folder(index + 5:endex - 1)];
mkdir(new_save_folder);
catch
    new_save_folder = ['D',new_save_folder(2:end)];
    mkdir(new_save_folder);
end
files = dir('*avg.mat');

notes = gather_notes(data_folder);

for i = 1:numel(files)
    load(files(i).name,'cdata');
    filename = [new_save_folder,files(i).name(1:end-3),'json'];
    save_cdata_as_json(filename, cdata, notes);
end
