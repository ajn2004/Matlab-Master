function im1 = color_it(i1,cmap,varargin)
if nargin == 3
    clims = varargin{1};
else
    clims = [min(i1(:)), max(i1(:))];
end
[m,n] = size(i1);
im1 = i1;
[mc,nc] = size(cmap);
si1 = (i1 - clims(1))/(clims(2)-clims(1));
si1 = si1.*(si1>0);
ind = find(si1 > 1);
si1(ind) = 1;
for i = 1:m
    for j = 1:n
        ind = round((mc-1)*si1(i,j))+1;
        cind = cmap(ind,:);
        for k = 1:nc
            im1(i,j,k) = cind(k);
        end
    end
end