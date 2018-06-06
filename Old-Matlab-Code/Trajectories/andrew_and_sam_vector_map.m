% function trajectory_map_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A a function to create a vector plot from the trajectories found in single
% molecule localization analysis
% AJN 7/12/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whitebg([0 0 0]);

clear i j
figure
hold on
xf_um=xf_all*q;
yf_um=yf_all*q;
xmin=0;
xmax=500;
ymin=0;
ymax=500;
time_per_frame=1/690;
cmap1=colormap(hot(10));        % Creates an RGB color matrix to color paths
cmap1(1,:)=[0.5 0.7 1];         % Defines Orange as fastest
Dbins=[0,.1,0.2,0.5,1,2,5,10,20,50,100];
nDbins=length(Dbins);

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
            for j=1:nDbins                      
                if Ds_all(i)<Dbins(j)
                    colorval=cmap1(j-1,:);
                    mobtype=j-1;
                    break
                end
            end
            if mobtype>1
                plot(xf_um(trajectories(i,1:nt)),yf_um(trajectories(i,1:nt)),'Color',colorval);
            else
                plot(xf_um(trajectories(i,1:nt)),yf_um(trajectories(i,1:nt)),'Color',colorval,'LineWidth',5);
            end
        end
%         waitbar(i/numel(trajectories(:,1)));
end
title('Trajectories');
xlabel('Position in um');
ylabel('Position in um');
hold off

figure
all_angs = 0;
%% Manipulation of Trajectories to measure angles
for i=1:numel(trajectories(:,1))    %cycles through each row of the trajectories
    nt=numel(trajectories(i,trajectories(i,:)>0)); %Cycles through all nonzero columns
    xset=xf_um(trajectories(i,1:nt));               % Creates a set of x variables
    yset=yf_um(trajectories(i,1:nt));               % Creates a set of y variables
    xdiff = diff(xset);
    ydiff = diff(yset);
    clear angle
    
    %this loop cycles through all the diff elements and corrects for
    %ambiguities in sign of angle that arises when using atan function
    for j = 1:numel(ydiff)
        if xdiff(j) >0 && ydiff(j) >0
           angle(j,1) = atan(ydiff(j)/xdiff(j));
        elseif xdiff(j) < 0 && ydiff(j) > 0
           angle(j,1) = pi+atan(ydiff(j)/xdiff(j));
        elseif xdiff(j) < 0 && ydiff(j) < 0
           angle(j,1) = pi + atan(ydiff(j)/xdiff(j));
        elseif xdiff(j) > 0 && ydiff(j) < 0
            angle(j,1) = atan(ydiff(j)/xdiff(j)) + 2*pi;
        end
    end
    all_angs = vertcat(all_angs,angle);
end

% define angle of cluster / region
param = polyfit(xf_um,yf_um,1);

% This value is determined for the cluster I am looking at, it assumes a
% preferential direction and should not be used if you do not have a
% preferential angle
% clus_ang = pi+atan(param(1));

% No implied preferntial direction of cluster
clus_ang = atan(param(1));

corr_angs = clus_ang-all_angs(2:numel(all_angs));
deg_ang = corr_angs*360/(2*pi);
hist(corr_angs.*360./(2*pi),-360:10:360);
xlabel('Angle in Degrees from best fit');
ylabel('Frequency in Number');
