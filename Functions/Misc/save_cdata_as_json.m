function save_cdata_as_json(file_name, cdata)
json_string = cdata2json(cdata);
fid = fopen(file_name, 'w');
fwrite(fid,json_string,'char');
fclose(fid);
end

