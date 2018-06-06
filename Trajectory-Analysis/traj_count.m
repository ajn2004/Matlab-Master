% traj_count
% a quick script to count the number of trajectories over different dmax
% files and display the resulting graph of traj. number vs dmax

files = dir('*.mat');
trajs = [];
lng = [];
for i = 1:numel(files)
    load(files(i).name,'trajec');
    trajs(i) = numel(trajec);
    ind1 = strfind(files(i).name,'dast') +5;
    ind2 = strfind(files(i).name,'nm') -1;
    lng(i) = str2num(files(i).name(ind1:ind2));
end

plot(lng,trajs,'.');