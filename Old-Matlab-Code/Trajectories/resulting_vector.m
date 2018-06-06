%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will manipulate data in the trajectory data set to give 1
% resulting vector displayed next to an arbitrary unit vector pointing in
% the x direction as measured by the image orientation. 
% AJN 10/3/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Loading the data set
[fname, fpath] = uigetfile('*traj.mat','Choose a trajectory file');
load(strcat(fpath,fname));

%% Declare Variables
ms_per_frame=1000/690;
xf_um = xf_all*q;
yf_um = yf_all*q;

%% Manipulation of Trajectories
for i=1:numel(trajectories(:,1))    %cycles through each row of the trajectories
    nt=numel(trajectories(i,trajectories(i,:)>0)); %Cycles through all nonzero columns
    xset=xf_um(trajectories(i,1:nt));               % Creates a set of x variables
    yset=yf_um(trajectories(i,1:nt));               % Creates a set of y variables
    xdiff = diff(xset);
    ydiff = diff(yset);
    sum_diff(i,1) = sum(xdiff);
    sum_diff(i,2) = sum(ydiff);
end


param = polyfit(xf_um,yf_um,1);
clus_ang = atan(param(1));

final(1,1) = sum(sum_diff(:,1));
final(1,2) = sum(sum_diff(:,2));
unit(1,1) = 1;
unit(1,2) = 0;
figure
hold on
quiver(0,0,final(1,1),final(1,2),'r');
quiver(0,0,-cos(clus_ang),-sin(clus_ang));
title('Net Displacement');
xlabel('Displacement in um');
ylabel('Displacement in um');
axis equal
hold off

% figure
% hold on
% quiver(0,0,final(1,1)/ms_per_frame,final(1,2)/ms_per_frame,'r');
% quiver(0,0,unit(1,1)/ms_per_frame,unit(1,2)/ms_per_frame,'b');
% title('Net Velocity');
% xlabel('velocity in um/ms');
% ylabel('velocity in um/ms');
% axis equal
% hold off