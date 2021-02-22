function notes = gather_notes(data_folder)
endex = strfind(data_folder,'\Analysis');
notes_path = data_folder(1:endex);
notes_file = [notes_path,'notes.txt'];

file = fopen(notes_file, 'r');
total_notes = fgetl(file);

while ischar(total_notes) % scan through line by line

total_notes = fgetl(file);
% look for red species info
found_str = strfind(total_notes,'red:');
if found_str ~= -1
    notes.red = total_notes(6:end);
end
found_str = strfind(total_notes,'orange:');
if found_str ~= -1
    notes.orange = total_notes(9:end);
end
found_str = strfind(total_notes,'info:');
if found_str ~= -1
    notes.info = total_notes(7:end);
end
end
fclose(file);


end