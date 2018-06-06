% Make Trajs
% A script to simulate data coming from da_storm with definite trajectories

xf_all = 1:0.01:3;
yf_all = xf_all;
framenum_all = 1:numel(xf_all);
q = 0.133;
dmax = 30;
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

