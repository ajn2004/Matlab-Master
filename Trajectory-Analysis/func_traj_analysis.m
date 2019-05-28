function [trajec] = func_traj_analysis(bname,dmax)
% Function traj analysis
% rewriting traj analysis to be called as a function
load([bname]); % loads the file selected
mkdir('traj');
trajec = struct('t',{[]}); % initialize trajectory variable
framenum_all = framenumber;
foll = zeros(numel(framenum_all),1);
xf_all = ncoords(:,1);
yf_all = ncoords(:,2);
zf_all = func_shift_correct(ncoords(:,3)*q,framenumber,1);
% zf_all = ncoords(:,3);
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
                    dist(k) = q*1000*((xf_all(cind(j)) - xf_all(fodex(k))).^2 + (yf_all(cind(j)) - yf_all(fodex(k))).^2).^0.5;
%                     dist(k) = q*1000*((xf_all(cind(j)) - xf_all(fodex(k))).^2 + (yf_all(cind(j)) - yf_all(fodex(k))).^2 + (zf_all(cind(j)) - zf_all(fodex(k))).^2).^0.5;  
                end
                disto = q*1000*((xf_all(cind(j)) - xf_all(cind)).^2 + (yf_all(cind(j)) - yf_all(cind)).^2).^0.5;
%                 disto = q*1000*((xf_all(cind(j)) - xf_all(cind)).^2 + (yf_all(cind(j)) - yf_all(cind)).^2 + (zf_all(cind(j)) - zf_all(cind)).^2).^0.5; 
%                 
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
for i = 1:numel(foll)
    traj = [];
    if foll(i) ~= 0
        traj = i;
        [traj, foll] = traj_connect(foll,traj, foll(i));
    trajec(count) = struct('t', traj);
  
    count = count+1;
    end
    
end
clear flag fodex i j k disto foll cind xf_all yf_all zf_all traj framenum_all
clear fname count
save(['traj\',bname(1:end-4),'_',num2str(dmax),'nm_traj.mat']);
end
