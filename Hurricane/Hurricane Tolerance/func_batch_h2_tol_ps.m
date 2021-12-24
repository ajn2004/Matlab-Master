function func_batch_h2_tol_ps(folder, file)
% batch routine for files
try
    load('C:\Users\ajnel\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
catch
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\Tolfile.mat');
end
fname = [folder, file];
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
color = 'red';
% sx_index = cal.red.z0s > tol.r.zlims(1) & cal.red.z0s < tol.r.zlims(2);
% sx_data = [cal.red.sx(sx_index);cal.red.sy(sx_index)].';
% nn_index = knnsearch(sx_data,[cdata.(color).fits(:,4),cdata.(color).fits(:,5)]);
% for i = 1:numel(cdata.(color).xf)
%     cdata.(color).distance(i,1) =  ((cdata.(color).fits(i,4) - sx_data(nn_index(i),1)).^2 + (cdata.(color).fits(i,5) - sx_data(nn_index(i),2)).^2)^0.5;
% end
%% defining statistical red
if channel_flag == 1 || channel_flag == 3 % load red tolerances
    cdata = apply_tolerance_to_cdata(cdata, tol, 'red');
end
%%  defining statistical orange
if channel_flag == 1 || channel_flag == 2 % load orange tolerances
    cdata = apply_tolerance_to_cdata(cdata, tol, 'orange');
end

% Scream here for sanity
save([folder(1:end-4),'Tol\',file(1:end-4),'_tol.mat']);
end