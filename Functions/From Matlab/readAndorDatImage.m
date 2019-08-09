function [andorImage, bFileDoesNotExist] = readAndorDatImage (rawDataFileName, param_file_name)
%param_file_name = '/Users/kannanuv/Downloads/metaSpool.txt';
%rawDataFileName = '/Users/kannanuv/Downloads/spool_X0589.dat';
%convertedImageFile = [rawDataFileName(1:end-4) '.tif'];

%% parse the paramFile here
param_ch1 = struct('fileName', param_file_name);
param_file_ptr = fopen (param_file_name, 'r');
tokeniser = ':';
line = fgetl (param_file_ptr);
while (~(line == -1))
    paramValue = line (1:(strfind (line, tokeniser) - 1));
    paramName = line ((strfind (line, tokeniser) + 1):end);
    if (strfind (paramName, 'AOI'))
        paramValue = str2double(paramValue);
    end
    if (strfind (paramName, 'Bytes'))
        paramValue = str2double(paramValue);
    end
    param_ch1 = setfield (param_ch1, paramName, paramValue);
    nToken = length (strfind (line, tokeniser)) + 1;
    line = fgetl (param_file_ptr);
end
fclose (param_file_ptr);

%% Call C function to readAndor dat file
filePtr = fopen (rawDataFileName);
if (filePtr == -1)
    bFileDoesNotExist = 1;
else
    bFileDoesNotExist = 0;
    fclose (filePtr);
end
if (bFileDoesNotExist)
    andorImage = zeros (param_ch1.AOIWidth, param_ch1.AOIHeight);
else
    andorImage = readAndorDatFile(rawDataFileName, param_ch1.AOIWidth, param_ch1.AOIHeight, param_ch1.AOIStride, param_ch1.ImageSizeBytes);
end
andorImage = uint16 (andorImage);
%imwrite (uint16(andorImage), convertedImageFile);