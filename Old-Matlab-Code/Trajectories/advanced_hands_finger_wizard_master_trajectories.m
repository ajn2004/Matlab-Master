%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a master batch file for testing trajectories.
% The script will perform the trajectory analysis at various dmax values
% then offer a plot to show the response of diffusion coefficient and
% number of trajectories found vs. dmax
% weight_ave gives out the average diffusion coefficient of a particular
% dmax value
% created 5/1/14 AJN
% This codes was built off of previously used codes

clear all
clc

dmax=0.3;   %max distance (in um) in next frame to be considered as a pair
dfactor=2;  %if next-nearest is within dmax*dfactor in the next frame, don't count as a pair
dfactor2=2; %if any within dmax*dfactor2 in the same frame, don't count as a pair
microns_per_pixel=0.133; %camera pixel size
exp_time = 0.02;    % frame exposure time in seconds
savename = 'auto_focus_test_dast.mat';

dirname='C:\Users\AJN Lab\Dropbox\Data\3-29-18 astig_cals\Analysis\';

%% Pair Identification
%dirname='C:\Users\Sam\Desktop\2Color_LiveCell\Analysis\05-14-2010\Analysis\CytoD_Dendra2-HA_PAmCherry-Actin\tolerances052410\';
cd(dirname);
files=dir(fullfile(dirname,'*_tol.mat'));
N_files=size(files);
mol_type=1;
for dmax = 0.01:0.01:1         % this for loop chooses various sizes of dmax to analyze

for cont1=1:N_files(1);

    this_file=files(cont1).name;
    this_file=this_file(1:end-4);
    matfile=strcat(dirname,this_file,'.mat');
%     index=find(this_file=='_tol');
    save_name=strcat(dirname,'',this_file,'_',num2str(dmax*1000),'nm_connect.mat');

    %data_all=[index_c index_n frame_c frame_n xc yc NN_dx NN_dy NN_dist_c NN_dist_n NNN_dist_n a0_c a0_n]
    [q,dist_data,xf_all,yf_all,framenum_all,a0_all,total_molecules,nrat,green_sum,red_sum]=next_frame_dist_func_TG1_2color2(matfile,microns_per_pixel,dmax,dfactor,dfactor2,mol_type);
 %load(matfile,'xf_all','yf_all','framenum_all','a0_all','total_molecules','nrat','green_sum','red_sum');
 %stop
    save(save_name,'q','xf_all','yf_all','nrat','red_sum','green_sum','dist_data','microns_per_pixel','dmax','dfactor','dfactor2');    
end
end

%% Trajectory creator
files_con=dir('*_connect.mat');
N_files_con=size(files_con);

for i=1:N_files_con(1);

    this_file=files_con(i).name;
    this_file=this_file(1:end-4);
    matfile=strcat(dirname,this_file,'.mat');
%     index=find(this_file=='_');
    save_name=strcat(matfile(1:end-4),'_traj.mat');

    trajectories_func;   
end

%% Trajectory analysis
File_Info=dir('*_traj.mat');
nfiles=length(File_Info);
Ds_all = zeros(N_files(1,1),round(5/0.05));
num_o_traj = zeros(N_files(1,1),round(5-0.05));
for j=1:nfiles
    clear ds_all_ave num_con dmax

    mat_name=File_Info(j).name;
    [ds_all_ave, num_con, dmax]=trajectory_analysis_func(dirname,mat_name,exp_time);
    if exist('dmax_ind','var')                           
        index = find(ismember(dmax_ind,dmax));
        if 1-isempty(index)
            iddqd = find(Ds_all(:,index)==0,1);
            Ds_all(iddqd,index) = ds_all_ave;
            num_o_traj(iddqd,index) = num_con;
        else
            dmax_ind = vertcat(dmax_ind,dmax);      
            Ds_all(1,numel(dmax_ind)) = ds_all_ave;
            num_o_traj(1,numel(dmax_ind)) = num_con;
        end
    else
        dmax_ind = dmax;
            Ds_all(1,1) = ds_all_ave;
            num_o_traj(1,1) = num_con;
    end
end

%% creates weighted averages of diffusion coefficients based on number of molecules measured in each file
for lk = 1:numel(dmax_ind)
    weight_ave(lk) = sum(num_o_traj(:,lk).*Ds_all(:,lk))/(sum(num_o_traj(:,lk)));
    tot_num(lk) = sum(num_o_traj(:,lk));
    err_diff(lk) = std(Ds_all(:,lk));
end

%% Representing the data
figure
subplot(2,1,1)
errorbar(dmax_ind,weight_ave,err_diff,'.b');
title('Diffusion Coefficient vs. dmax');
xlabel('dmax in um');
ylabel('Diffusion Coefficient in um^2/s');

subplot(2,1,2)
plot(dmax_ind,tot_num,'.b');
title('Total Number of Trajectories vs. dmax');
xlabel('dmax in um');
ylabel('Number of trajectories found');

%% Makes an estimation of the diffusion coefficient based on total number of trajectories found
INdex = logical(1-isnan(weight_ave));
diffuse_guess = sum(weight_ave(INdex).*tot_num(INdex))/sum(tot_num(INdex));
save(savename,'weight_ave' , 'tot_num', 'num_o_traj', 'Ds_all','dmax_ind', 'err_diff');