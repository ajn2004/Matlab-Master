%quick_cor

for i = 1:numel(tinfo)
    load(tinfo(i).name);
    fluor2(:,i) =fluor(:);
end

aves = mean(fluor2,1);

Es = mean((fluor2(:,1)-aves(1)).*(fluor2(:,2)-aves(2)));
corr = Es/(std(fluor2(:,1))*std(fluor2(:,2)));