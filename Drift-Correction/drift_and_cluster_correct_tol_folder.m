function drift_and_cluster_correct_tol_folder(folder,color, framechunk)
files = dir([folder, '*tol*']);
for i = 1:numel(files)
    try
    load([files(i).folder,'\',files(i).name]);
    cdata = clean_cdata(cdata);
    cdata = model_drift_correction(cdata,color,framechunk);
%     cdata = cluster_clean(cdata, 0.15, 5);
    save([folder(1:end-4),'Fin\',files(i).name(1:end-7), 'fin.mat'],'cdata','cal');
%     delete(files(i).name);
    catch lsterr
        disp(lsterr)

    end
end