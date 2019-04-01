% Quick Look
%
% This will quickly organize the trajectory data to show number of
% trajectories as a function of dmax

clearvars; close all; clc; % cleanup
files = dir('*traj.mat'); % grab all file info

% Preallocation
trajs = zeros(numel(files),1);
trajl = [];
dms = trajs;
for i = 1:numel(files) % looping over files
    load(files(i).name,'trajec','q','ncoords','dmax','llv'); % load file
    id = strfind(files(i).name,'pre'); % Find index if 
    if isempty(id) % determining pre/post 
        indy(i) = 0;
    else
        indy(i) = 1;
    end
    trajs(i) = numel(trajec); % Record the number of of trajectories found
    trajl(i) = 0; % Preallocate trajectory length
    clear dista % reset distance variable
    dx =[];
        dy =[];
        dz =[];
    for j = 1:trajs(i) % loop over trajects
        trajl(i) = trajl(i) + numel(trajec(j).t); % Record trajectory length
        ind = trajec(j).t; % Grab trajectories
        
        for k = 1:(numel(ind)-1) % Loop over trajectory length
            dista(j) = 0; % Reseting dista variable
            for l = 1:3 % Lazily loop over coords
                dista(j) = dista(j) + (ncoords(ind(k),l) - ncoords(ind(k +1),l))^2;
            end
            dista(j) = q*dista(j)^0.5; % adding in squares needs to square root it later
            dx = [dx;q*(ncoords(ind,1)-mean(ncoords(ind,1)))];
            dy = [dy;q*(ncoords(ind,2)-mean(ncoords(ind,2)))];
            dz = [dz;q*(ncoords(ind,3)-mean(ncoords(ind,3)))];
        end
        
    end
    msx(i) = (sum((dx).^2)/(numel(dx)-1)).^0.5; % Standard deviation of steps
    msy(i) = (sum((dy).^2)/(numel(dx)-1)).^0.5;
    msz(i) = (sum((dz).^2)/(numel(dx)-1)).^0.5;
    mdist(i) = mean(dista(dista>0)); % Record mean distance of a trajectory step
    trajl(i) = trajl(i)/trajs(i); % Record trajcetory length / number of trajectories
    dms(i) = dmax; % record maximum distance
    trajp(i) = trajs(i)./numel(llv); % record trajectories / number of molecules
    clear trajec dmax
end

plot(dms,trajs,'.');
title('Number of trajectories')
xlabel('Distance cut off in nm');
ylabel('Number of Trajectories');
figure
plot(dms,trajl,'.');
title('Total Trajectory Length normalized to total trajectories')
xlabel('Distance Cut off in nm')
ylabel('Length/Trajs')
figure
plot(dms,mdist,'.');
title('Average Trajectory Step')
xlabel('Distance Cut Off in nm')
ylabel('Average Step')
figure
indy = logical(indy);
plot(dms(indy),(msx(indy).^2+msy(indy).^2).^0.5*1000,'.r');
hold on
plot(dms(logical(1-indy)),(msx(logical(1-indy)).^2+msy(logical(1-indy)).^2).^0.5*1000,'.b');
% plot(dms,msy*1000,'.');
title('Lateral \sigma step')
xlabel('Distance cut off nm')
legend('Prestimulation','Poststimulation')
ylabel('Lateral \sigma in nm')
figure
plot(dms,msz*1000,'.')
title('Axial \sigma step')
xlabel('Distance cut off nm')
ylabel('Axial \sigma in nm')