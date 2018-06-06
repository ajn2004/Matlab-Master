function cmdstrng  = quicksave(savename1,extension, names,varys)
% this is a quick script that allows the user to list names over variables
% to save in the 'varys' variable, and the names of other variables to be
% saved and mentioned in the filename in the variable 'names'
cmdstrng = ['save(','''',savename1];
varstring = ');';
for i = 1:numel(varys)
    varstring = [',','''',varys{i},'''',varstring];
end

for i = 1:numel(names)
    cmdstrng = [ cmdstrng,'_',names{i}];
    varstring = [',','''',names{i}, '''', varstring];
end
    cmdstrng = [cmdstrng,extension,'''',varstring];
%     eval(cmdstrng)
end