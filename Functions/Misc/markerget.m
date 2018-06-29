function output = markerget(str,marker)
output = [];
ind1 = strfind(str,marker);
ind2 = -1 + ind1 + strfind(str(ind1:end),'_');
if isempty(ind2) == 1
    ind2 = numel(str)+1;
end
output = str(ind1+1:ind2(1)-1);
output = str2num(output);

