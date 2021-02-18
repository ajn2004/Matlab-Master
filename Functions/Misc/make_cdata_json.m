clearvars
mkdir('jsons')
files = dir('*avg.mat');

for i = 1:numel(files)
    load(files(i).name,'cdata');
    filename = ['jsons\',files(i).name(1:end-3),'json'];
    save_cdata_as_json(filename, cdata);
end
