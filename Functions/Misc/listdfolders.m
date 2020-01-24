function F = listdfolders(path)
% This function will return an object F listing all folders under path
f = genpath(pwd);


ind = strfind(f(1:end),';');
for i = 1:numel(ind)
    if ~exist('F')
        F{1} = [f(i:ind(i)-1),'\'];
    else
        F{numel(F)+1} = [f(ind(i-1)+1:ind(i)-1),'\'];
    end
end
