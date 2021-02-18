clearvars
base_folder = 'G:\Dropbox\Andrew & Tim shared folder\Analyzed Data\';
data_folder = pwd;
index = strfind(data_folder, 'Data');
endex = strfind(data_folder, 'Analysis');

new_save_folder = [base_folder, data_folder(index + 5:endex - 1)];
mkdir(new_save_folder);
files = dir('*avg.mat');

for i = 1:numel(files)
    load(files(i).name,'cdata');
    filename = [new_save_folder,files(i).name(1:end-3),'json'];
    save_cdata_as_json(filename, cdata);
end
