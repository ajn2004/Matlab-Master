% Reading the Serial Port and Displaying the output
scom = 'Com3';

% Declare serial object
s = serial(scom);
% Open it
fopen(s);

% Doing stuff here
cnt = 1;
while true
    try
        tic
    str = fscanf(s);
    ind = strfind(str,',');
    v1(cnt) = str2num(str(1:ind(1)-1));
    v2(cnt) = str2num(str(ind(1)+1:ind(2)-1));
    count(cnt) = str2num(str(ind(2)+1:end));
    cnt = cnt +1;
    toc
    catch lsterr
    end
end

% Close the serial object
fclose(s);