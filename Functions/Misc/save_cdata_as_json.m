function save_cdata_as_json(file_name, cdata, notes)
json_string = cdata2json(cdata, notes);
fid = fopen(file_name, 'w');
fwrite(fid,json_string,'char');
fclose(fid);
end

