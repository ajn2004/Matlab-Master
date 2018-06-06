% Example start up code for matlab utilization
% This will have matlab add all folders and subfolders associated with the
% gitHub downloaded files and allow it to access any scripts associated
% with my codes

ghpath = '[File Path to GitHub Folder]';
addpath(ghpath);
addpath(genpath(ghpath));
