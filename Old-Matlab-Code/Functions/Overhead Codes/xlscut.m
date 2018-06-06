% XLS Convertion for Ryan
function s = xlscut(fname)
s =[];
[finfo, sheets] = xlsfinfo(fname);

for i = 1:numel(sheets)
    str = ['s.',sheets{i}, ' = xlsread(fname,i);'];
    eval(str);
end