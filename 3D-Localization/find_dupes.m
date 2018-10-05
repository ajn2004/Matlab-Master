function ind = find_dupes(cents,olap)
% Find duplicates that overlap by too much
x = cents(:,1);
y = cents(:,2);
ind = [];
for i = 1:numel(x)
    d2 = (x-x(i)).^2 + (y-y(i)).^2;
    d = sort(d2.^0.5);
    if d(2) < olap
        ind = [ind;i];
    end
end
    