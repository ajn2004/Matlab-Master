% This is a function to analyze the data created in our trajectory analysis
% program
% Created 5/1/14 AJN

function[Ds_all_ave, num_con, dmax] = trajectory_analysis_func(fpath,fname,time_per_frame)
load([fpath,'\',fname]);
clear i j
xf_um=xf_all*q;
yf_um=yf_all*q;
xmin=0;
xmax=500;
ymin=0;
ymax=500;
if exist('trajectories','var')
for i=1:numel(trajectories(:,1))    %cycles through each row of the trajectories
    nt=numel(trajectories(i,trajectories(i,:)>0)); %Cycles through all nonzero columns
    xset=xf_um(trajectories(i,1:nt));
    yset=yf_um(trajectories(i,1:nt));
        if xset(1)>xmin && xset(1)<xmax && yset(1)>ymin && yset(1)<ymax
            xdisp=diff(xset);                   % Creates a vector of differences in x values
            ydisp=diff(yset);                   % Creates a vector of differences in y values
            sd=xdisp.*xdisp+ydisp.*ydisp;       % Creates a vector of distances between localizations
            msd=mean(sd);                       % Average distance traveled
            Ds_all(i)=msd/(4*time_per_frame);   % Diffusion Coefficient of particle

        end
%         waitbar(i/numel(trajectories(:,1)));
end
Ds_all_ave = mean(Ds_all);
num_con = numel(trajectories);
else
    Ds_all_ave = 0;
    num_con = 0;
end
