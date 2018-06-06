% Analyze the results
% startitup;

[fname, fpath] = uigetfile('*mat');

cd(fpath);
finfo = dir('*mat');

for i = 1:numel(finfo)
    load(finfo(i).name);
    nums = numel(xf_all);
    for j = 1:nums
        ind = find(number == j);
        xcm = sum(xf_all(ind))/numel(ind);
        ycm = sum(yf_all(ind))/numel(ind);
        dists = ((xf_all(ind) - xcm).^2 + (yf_all(ind) - ycm).^2).^0.5;
        ddists = (diff(xf_all(ind,1)).^2 + diff(yf_all(ind)).^2).^0.5;
        aveddists(j,i) = mean(ddists)*156;
        avedist(j,i) = mean(dists)*156;
        clear dists
    end
end