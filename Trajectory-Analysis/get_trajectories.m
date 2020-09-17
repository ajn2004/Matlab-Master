function [trajec] = get_trajectories(data,frames,dmax)
% Function traj analysis
% rewriting traj analysis to be called as a function

% trajec = struct('t',{[]}); % initialize trajectory variable
framenum_all = frames;
foll = zeros(numel(framenum_all),1);
xf_all = data(:,1);
yf_all = data(:,2);

% loop over all frames to build the connections
for i = 1:max(framenum_all)
    clear dist
    % i is the framenumber, so we want to look at
    if i ~= max(framenum_all) % special conditions for frame 1 because there is no previous
        cind = find(framenum_all == i); % current index for loop i
        fodex = find(framenum_all == i+1); % index of all molecules on following frame
        if ~isempty(cind) && ~isempty(fodex)
            for j = 1:numel(cind)
                clear dist disto
                for k = 1:numel(fodex)
                    try
                        dist(k) = 0;
                        for l = 1:numel(data(1,:))
                        dist(k) = dist(k) +(data(cind(j),l) - data(fodex(k),l)).^2;
                        end
                        dist(k) = dist(k)^0.5;
                    catch
                        trajec = 'new_baby';
                    end
                end

%                     disto = ((xf_all(cind(j)) - xf_all(cind)).^2 + (yf_all(cind(j)) - yf_all(cind)).^2).^0.5;
                    disto = 0;
                    for l = 1:numel(data(1,:))
                        disto = disto +(data(cind(j),l) - data(cind,l)).^2;
                    end
                    disto = disto.^0.5;
                flag = find(disto > 0 & disto <= 2*dmax);
                if min(dist) <= dmax && sum(dist <= 2*dmax) < 2 && isempty(flag)
                    foll(cind(j)) = fodex(find(dist == min(dist)));
                end
            end
        end
    end
end

%% loop over conections to build trajectory variable
count = 0;
trajec = xf_all * 0;
for i = 1:numel(foll)
    if foll(i) ~= 0 && trajec(i) == 0
        count = count + 1;
        trajec = label_follows(foll,trajec,count,i);
    end
    
end
end

function trajec = label_follows(foll,trajec,trajectory_number,current)
    trajec(current) = trajectory_number;
    if foll(current) ~= 0     
        trajec = label_follows(foll,trajec,trajectory_number,foll(current));
    end
end
