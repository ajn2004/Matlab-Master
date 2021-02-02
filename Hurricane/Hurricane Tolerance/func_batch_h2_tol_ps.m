function func_batch_h2_tol_ps(fname)
% batch routine for files
try
load('C:\Users\ajnel\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
catch
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
end
load(fname,'cdata','cal','pixw','q');

% Identify number of colors we're working with
flag = 0;
try
    xf = cdata.red.xf; % if red loads successfully, start red tolerances
    flag = 1; % if red exists assume 2-color, if orange also exists flag == 1
    
catch lsterr
end

try
    xf = cdata.orange.xf; % if orange loads successful, start orange tolerances
    if flag == 0 % if no red data this is an 'orange only'
        flag = 2;
    end
    
catch
    flag = 3; % Orange channel failed indicates must have red data, or will cause error
end
channel_flag = flag;
% channel_flag
% if channel flag = 1, dual color
% if channel flag = 2, orange only
% if channel flag = 3, red only
tol.q = cal.q;
%% defining statistical red
if channel_flag == 1 || channel_flag == 3 % load red tolerances
    cdata = apply_tolerance_to_cdata(cdata, tol, 'red');
end
%%  defining statistical orange
if channel_flag == 1 || channel_flag == 2 % load orange tolerances
    cdata = apply_tolerance_to_cdata(cdata, tol, 'orange');
end

% Scream here for sanity
save([fname(1:end-4),'_tol.mat']);
end