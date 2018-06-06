clear all
clc

dirname='I:\Data\4-18-14 2color ffpalm\Analysis\30 ms\Toleranced\';
%dirname='C:\Users\Sam\Desktop\2Color_LiveCell\Analysis\05-14-2010\Analysis\CytoD_Dendra2-HA_PAmCherry-Actin\tolerances052410\';
cd(dirname)
files=dir('*_connect.mat');
N_files=size(files);

for i=1:N_files(1);
    this_file=files(i).name;
    this_file=this_file(1:end-4)
    matfile=strcat(dirname,this_file,'.mat');
%     index=find(this_file=='_');
    save_name=strcat(matfile(1:end-4),'_traj.mat');

    trajectories_func;   
end