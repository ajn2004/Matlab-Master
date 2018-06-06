files = dir('*dast*');

for i = 1:numel(files)
    load(files(i).name);
    zf_nm = getdz(sigx_all,sigy_all);
    save(files(i).name);
end