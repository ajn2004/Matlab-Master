function dmin = get_mins(f)
[m,n] = size(f);
mins(1) = f(1);
inds = 1;
for i = 2:max([m,n])
    if f(i) <= mins(end)
        mins(i) = f(i);
        inds = [inds;i];
    end
end
dmin = spline(inds,mins,1:max([m,n]));