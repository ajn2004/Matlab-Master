clear all
clc

dmax=0.3;   %max distance (in um) in next frame to be considered as a pair
dfactor=2;  %if next-nearest is within dmax*dfactor in the next frame, don't count as a pair
dfactor2=2; %if any within dmax*dfactor2 in the same frame, don't count as a pair
microns_per_pixel=0.106; %camera pixel size

dirname='I:\Data\4-18-14 2color ffpalm\Analysis\30 ms\Toleranced\';
%dirname='C:\Users\Sam\Desktop\2Color_LiveCell\Analysis\05-14-2010\Analysis\CytoD_Dendra2-HA_PAmCherry-Actin\tolerances052410\';

files=dir(fullfile(dirname,'*_tol.mat'));
N_files=size(files);
mol_type=1;
for dmax = 0.05:0.05:0.5
for i=1:N_files(1);
    this_file=files(i).name;
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

files_con=dir('*_connect.mat');
N_files_con=size(files_con);

for i=1:N_files_con(1);
    this_file=files_con(i).name;
    this_file=this_file(1:end-4)
    matfile=strcat(dirname,this_file,'.mat');
%     index=find(this_file=='_');
    save_name=strcat(matfile(1:end-4),'_traj.mat');

    trajectories_func;   
end
