% Wave Trajectories
% Display trajectories of wavelet analysis localizations
% close all
% xf = fits(:,1) + cents(:,1);
% yf = fits(:,2) + cents(:,2);
% N = fits(:,3);

dmax = 500;
% ind  = fits(:,3) > 70;
% ind  = N > 200;

xf_all = dcd(ind,1);
yf_all = dcd(ind,2);
zf_all = dcd(ind,3);
framenum_all = fms(fnum(ind)); % convert to absolute framenumber
trajec = struct('t',{[]}); % initialize trajectory variable
foll = zeros(numel(framenum_all),1);
% xc = fits(ind,1);
% yc = fits(ind,2);
% sdi1 = iloc(:,:,ind);
% iloc = sdi1;
%% This has been copied from 'traj_analysis.m'
% loop over all frames to build the connections
for i = 1:max(framenum_all)
    clear dist
    % This is a modified trajectory analysis. We know the frame the 'first'
    % molecule will appear on, and want to only look at subsquent ones
    if ismember(i,stimdex) % special conditions for stimulus frames because we don't want to connect w/ a previous frame
        cind = find(framenum_all == i); % current index for loop i
        fodex = find(framenum_all == i+1); % index of all molecules on following stim frames
        for j = 1:fps - bk_fms - 1
            fodex = [fodex,find(framenum_all == i +j+1)];
        end
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
count = 1;
for i = 1:numel(foll) % foll is an array of follows
    traj = []; % initiallize traj number
    if foll(i) ~= 0 % if there is a following index
        traj = i; % Start traj w/ current index
        [traj, foll] = traj_connect(foll,traj, foll(i)); % send all follows, the traj you're building, and the current follower
    trajec(count) = struct('t', traj);
  
    count = count+1;
    end
    
end

set_scale(std(dip1,1,3),q,4);
colormap('gray');
hold on
[~, ind1,~] = intersect(framenum_all,stimdex);
[~, ind2] = setdiff(framenum_all,stimdex);

zf_all = zf_all - min(zf_all);
plot3(xf_all(ind1),yf_all(ind1),zf_all(ind1),'.g')
plot3(xf_all(ind2),yf_all(ind2),zf_all(ind2),'.r')
for i = 1:numel(trajec)    
    inds = trajec(i).t;
    d(i) = q*((xf_all(trajec(i).t(1)) - xf_all(trajec(i).t(2))).^2 + (yf_all(trajec(i).t(1)) - yf_all(trajec(i).t(2))).^2 + (zf_all(trajec(i).t(1)) - zf_all(trajec(i).t(2))).^2).^0.5;
    plot3(xf_all(inds),yf_all(inds),zf_all(inds),'b');
end
legend('Releases','Post-stims','Trajectory')
hold off
axis image
figure
plot3(xf_all(ind1)*q,q*yf_all(ind1),q*zf_all(ind1),'g.')
hold on
plot3(xf_all(ind2)*q,q*yf_all(ind2),q*zf_all(ind2),'r.')
hold off
axis equal
xlabel('X-axis nm')
ylabel('Y-axis nm')
zlabel('Z-axisl nm')
d = d.';
