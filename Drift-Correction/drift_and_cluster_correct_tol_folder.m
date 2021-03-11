function drift_and_cluster_correct_tol_folder()
files = dir("*tol.mat");
for i = 1:numel(files)
    load(files(i).name);
    cdata = clean_cdata(cdata);
    cdata = model_drift_correction(cdata,'red',1000);
    cdata = cluster_clean(cdata, 0.15, 5);
    save([files(i).name(1:end-7), 'fin.mat'],'cdata','cal');
%     delete(files(i).name);
end