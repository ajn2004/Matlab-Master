% A script to batch trajectory analysis
clear all
clc

mat_dir = 'I:\Data\4-18-14 2color ffpalm\Analysis\30 ms\Toleranced\';


cd(mat_dir);
File_Info=dir('*_traj.mat');
nfiles=length(File_Info);

for i=1:nfiles
    clear ds_all_ave num_con dmax
    i
    mat_name=File_Info(i).name;
    [ds_all_ave, num_con, dmax]=trajectory_analysis_func(mat_dir,mat_name);
    if exist('dmax_ind','var')                           
        index = find(ismember(dmax_ind,dmax));
        if 1-isempty(index)            
            Ds_all(i,index) = ds_all_ave;
            num_o_mol(i,index) = num_con;
        else
            dmax_ind = vertcat(dmax_ind,dmax);
            Ds_all(i,numel(dmax_ind)) = ds_all_ave;
            num_o_mol(i,numel(dmax_ind)) = num_con;
        end
    else
        dmax_ind = dmax;
            Ds_all(i,numel(dmax_ind)) = ds_all_ave;
            num_o_mol(i,numel(dmax_ind)) = num_con;
    end
end

for lk = 1:numel(dmax_ind)
    weight_ave(lk) = sum(num_o_mol(:,lk).*Ds_all(:,lk))/(sum(num_o_mol(:,lk)));
    tot_num = sum(num_o_mol(:,lk));
end

plot(dmax_ind,weight_ave,'.b')
title('Diffusion Coefficient vs. dmax');
xlabel('dmax in um');
ylabel('Diffusion Coefficient in um^2/s');