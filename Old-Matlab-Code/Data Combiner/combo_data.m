% Combo Data
% Combines large data sets into a single data set to be used on toleranced
% files
%
% 2/7/18 AJN
try
while true
clear all
close all
clc

nfiles = 1;
mkdir('tot');

% Use a try catch loop to load as many files as necessary, handles variable
% number of files
try
    while true
        
        [fname_temp, fpath_temp] = uigetfile('*mat','Select a file in chronological order, cancel to continue');
        cd(fpath_temp);
        fname{nfiles} = cellstr(fname_temp);
        fpath{nfiles} = cellstr(fpath_temp);
        nfiles = nfiles +1;
    end
catch lasterr
    nfiles=nfiles - 1;
    disp(['Combining ', num2str(nfiles), ' files']);
end

% preallocate initial arrays

N_temp =[];
N_crlb_temp =[];
framenum_all_temp = [];
off_temp = [];
off_crlb_temp = [];
% lp_temp =[];
% lp2_temp =[];
xf_all_temp =[];
xf_crlb_temp =[];
yf_all_temp =[];
yf_crlb_temp =[];
sigx_all_temp = [];
sigx_crlb_temp = [];
sigy_all_temp = [];
sigy_crlb_temp = [];

frames_last_file(1) = 0;
llv_temp = [];
y_temp = [];
iloc_temp = [];
imfiles = [];
total_mol = 0;
%loop over all files combining relevant localization data
for jk = 1:nfiles
    load([char(fpath{jk}), char(fname{jk})]);
    if jk == 1
        tot_name = char(fname{jk}); 
    end
    frames_last_file(jk+1) = framenum_all(end); %keep track of number of frames in this data set
                   N_temp = [N_temp; N];
              N_crlb_temp = [N_crlb_temp; N_crlb];
                 off_temp = [off_temp; off_all];
            off_crlb_temp = [off_crlb_temp; off_crlb];
        framenum_all_temp = [framenum_all_temp; framenum_all + sum(frames_last_file(1:jk))];
                 llv_temp = [llv_temp;llv];
                   y_temp = [y_temp;y];
                iloc_temp = [iloc_temp,iloc];
              xf_all_temp = [xf_all_temp; xf_all];
              yf_all_temp = [yf_all_temp; yf_all];
             xf_crlb_temp = [xf_crlb_temp; xf_crlb];
             yf_crlb_temp = [yf_crlb_temp; yf_crlb];
            sigx_all_temp = [sigx_all_temp; sigx_all];
           sigx_crlb_temp = [sigx_crlb_temp; sigx_crlb];
            sigy_all_temp = [sigy_all_temp;sigy_all];
           sigy_crlb_temp = [sigy_crlb_temp;sigy_crlb];
                total_mol = total_mol + numel(yf_crlb);
              imfiles{jk} = char(fname{jk});
    clear N N_crlb framenum_all lp lp2 y zf_all iloc sigx_all sigx_crlb
    clear xf_all xf_crlb yf_all yf_crlb off_all off_crlb sigy_all sigy_crlb
    clear total_molecules imagefile imagfo
end

% assign temp variables to final arrays
N = N_temp;
N_crlb = N_crlb_temp;
off_all = off_temp;
off_crlb = off_crlb_temp;
framenum_all = framenum_all_temp;
xf_all = xf_all_temp;
xf_crlb = xf_crlb_temp;
sigx_all = sigx_all_temp;
sigy_all = sigy_all_temp;
sigx_crlb = sigx_crlb_temp;
sigy_crlb = sigy_crlb_temp;
zf_all = getdz(sigx_all,sigy_all)/q;
llv = llv_temp;
yf_all = yf_all_temp;
yf_crlb = yf_crlb_temp;
base_name = tot_name;
total_molecues = total_mol;
iloc = iloc_temp;
y = y_temp;
% clear the mess
clear N_temp off_all_temp off_crlb_all_temp framenum_all_temp iloc_temp y_temp zf_all_temp 
clear lp_temp lp2_temp xf_all_temp xf_crlb_temp yf_crlb_temp N_crlb_temp zf_crlb_temp
clear yf_all_temp sigx_all_temp sigx_crlb_temp sigy_all_temp sigy_crlb_temp
clear off_crlb_temp orig_N n_end n_fail_a0 n_fail_outbox n_start
clear plot_histogram save_mat_file savename mat_name mat_file mat_dir
clear imagfo display_image dat_dir add_dir an_dir
save([char(fpath{jk}),'\tot\',base_name(1:end-4),'_tot.mat'])
end
catch lsterr
end