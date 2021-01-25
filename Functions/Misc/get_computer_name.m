function path_name = get_computer_name()
try
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat');
    path_name = 'C:\Users\AJN Lab';
catch 
    load('C:\Users\ajnel\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat');
    path_name = 'C:\Users\ajnel';
end
end