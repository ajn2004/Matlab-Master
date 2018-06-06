%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpu trajectories
%
% a trajectory code to be run on cpu
%
% ajn 2/6/16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% User Variables
dmax = 100; % maximum distance in nm

%% Preliminary Jazz
[fname, fpath] = uigetfile('*.mat'); % allows user to grab file
load([fpath,fname]); % loads the file selected

trajec = struct('t',{[]}); % initialize trajectory variable
foll = zeros(numel(framenum_all),1);

% loop over all frames to build the connections
for i = 1:max(framenum_all)
    clear dist
    % i is the framenumber, so we want to look at
    if i ~= max(framenum_all) % special conditions for frame 1 because there is no previous
        cind = find(framenum_all == i); % current index for loop i
        fodex = find(framenum_all == i+1); % index of all molecules on following frame
        if ~isempty(cind) && ~isempty(fodex)
            for j = 1:numel(cind)
                for k = 1:numel(fodex)
                    dist(k) = q*1000*((xf_all(cind(j)) - xf_all(fodex(k))).^2 + (yf_all(cind(j)) - yf_all(fodex(k))).^2).^0.5;
                end
                disto = q*1000*((xf_all(cind(j)) - xf_all(cind)).^2 + (yf_all(cind(j)) - yf_all(cind)).^2).^0.5; 
                flag = find(disto > 0 & disto <= 2*dmax);
                
                if min(dist) <= dmax && sum(dist <= 2*dmax) < 2 && isempty(flag)
                    foll(cind(j)) = fodex(find(dist == min(dist)));
                end
            end
        end
    end
end

%% loop over conections to build trajectory variable
count = 1
for i = 1:numel(foll)
    traj = [];
    if foll(i) ~= 0
        traj = i;
        [traj, foll] = traj_connect(foll,traj, foll(i));
    trajec(count) = struct('t', traj);
  
    count = count+1;
    end
    
end
save([fname(1:end-4),'_traj.mat']);
%% Represent trajectories
figure
hold on
for i = 1:numel(trajec)
    plot3(xf_all(trajec(i).t), yf_all(trajec(i).t),zf_nm(trajec(i).t));
end
hold off