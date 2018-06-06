
close all
clear all
clc

[fname, fpath] = uigetfile('*.mat');
dist = 1;
cd(fpath);

finfo = dir('*t_neuro.mat');
truth_name = 'mock_data_pN_150_pS_1_bkg_4_mpf_5noverlap.mat';
N_fa = [];
N_mi = [];
off_p = [];
N_p =[];
off_mi = [];
off_fa = [];
sigs_mi = [];

for i = 1:numel(finfo)
    i
   [missed_mol(i,1), false_p(i,1), victory(i,1), tot_num(i,1), N_fat, N_mit,N_pt, off_pt, off_mit, off_fat, sigs_mit  ] = func_data_compare(finfo(i).name,truth_name, dist); 
   N_fa = [N_fa;N_fat];
   N_mi = [N_mi;N_mit];
   off_p = [off_p; off_pt];
   N_p =[N_p;N_pt];
   off_mi = [off_mi;off_mit];
   off_fa = [off_fa;off_fat];
   sigs_mi = [sigs_mi;sigs_mit];
   inde = findstr(finfo(i).name, 'it_');
   index(i,1) = str2num(finfo(i).name(inde+3:end-12));
   

end


plot(index,false_p./(false_p + victory), '.r', 'MarkerSize', 30)
hold on
plot(index,victory./tot_num, '.b',  'MarkerSize', 30)
legend('Perecentage of false positives');
legend('Perecentage of false positives','Percentage of molecules discovered');
xlabel('Iteration Number');
ylabel('%');
title('Percentages as a function of iteration number')