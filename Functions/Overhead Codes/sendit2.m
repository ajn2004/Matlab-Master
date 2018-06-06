function sendit2(fname)
if ~exist(fname) % If a correction folder doesn't already exist
    mkdir(fname) % make one
end