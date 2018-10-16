function ind = find_fm_dupes(cents,fm,olap)
% Find duplicates that overlap by too much
x1 = cents(:,1);
y1 = cents(:,2);
ind = [];
for i = 1:numel(x1)
    id = fm == fm(i);
    d2 = (x1(id)-x1(i)).^2 + (y1(id)-y1(i)).^2;
    d = sort(d2.^0.5);
    if numel(d) >1
    if d(2) < olap
        ind = [ind;i];
    end
    end
end
    