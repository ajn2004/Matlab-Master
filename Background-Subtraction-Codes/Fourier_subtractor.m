[fname, fpath] = uigetfile('*.tif');
cd(fpath);
finfo = dir('*.tif');

for i = 1:numel(finfo)
    disp(['Loading number ', num2str(i)])
    i1 = readtiff(finfo(i).name);
    disp(['Subbing number ', num2str(i)])

        i2= fsub(i1,10);

    disp(['Writing number ', num2str(i)])
    writetiff(i2, [finfo(i).name(1:end-4),'_bg.tif']);
    clear i2 i1
    clc
end