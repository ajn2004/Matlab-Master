% Batch program to merge duplicates after scan correction
clearvars

folder = {'D:\Dropbox\Data\9-23-20 fixed hek gpi-halo\Analysis\toleranced\scan_corrected'};
cd(folder{1});
mkdir('converged')
files = dir('*sc.mat');
for i  = 1:numel(files)
    load(files(i).name);
    try
        cdata.orange = converge_duplicates(cdata.orange,2, 1);
    catch
    end
    try
        cdata.red = converge_duplicates(cdata.red,2, 0);
    catch
    end
    save(['converged\',files(i).name(1:end-4),'_cv.mat'],'cdata','tol','cal');
end

