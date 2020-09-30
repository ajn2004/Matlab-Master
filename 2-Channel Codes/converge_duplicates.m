function cdata_out = converge_duplicates(cdata_in,cluster_number, color)
% We're going to assume that incoming cdata is post-color in that cdata_in = cdata.orange
% preallocate outgoing cdata
ind_r = abs(cdata_in.zf) >= 0.5 ; % Requires pre-corrected z-values
%remove out of bounds data
fn = fieldnames(cdata_in);
for i = 1:numel(fn)
    if strcmp(fn{i},'fits') ||strcmp(fn{i},'crlbs')
        cdata_in.(fn{i})(ind_r,:) = [];
    else
        cdata_in.(fn{i})(ind_r) = [];
    end
end
% perform trajectory analysis
%     color = {'red','orange'};
% color_index = 1;
xf = cdata_in.xf;
yf = cdata_in.yf;
zf = cdata_in.zf_smoothed;
frames = cdata_in.framenumber;
if color == 1 % 1 for orange 0 for red
    frames = (frames+1)/2; % orange frame correction
else
    frames = (frames )/2; % Red frame correction
end
trajectories = get_trajectories([xf,yf,zf],frames,0.07);

x_tot = [];
y_tot = [];
z_tot = [];
for m = 1:max(trajectories)
    ind = find(trajectories == m);
    x_tot = [x_tot; xf(ind) - mean(xf(ind))];
    y_tot = [y_tot; yf(ind) - mean(yf(ind))];
    z_tot = [z_tot; zf(ind) - mean(zf(ind))];
    
end


cluster_linkage_distance = std(x_tot) + std(y_tot); % We want 2x the avg std, which is 2*(stdx + stdy )/2 = what we have

% Cluster analysis and merging
clusters = dbscan([xf,yf],cluster_linkage_distance,cluster_number);
% Preallocate array ,
ind = clusters == -1; % -1 are locs that are 'unclustered'
cdata_out = [];
for i = 1:numel(fn)
    if strcmp(fn{i},'fits') ||strcmp(fn{i},'crlbs')
        cdata_out.(fn{i})(ind,:) = cdata_in.(fn{i})(ind,:);
    else
        cdata_out.(fn{i})(ind) = cdata_in.(fn{i})(ind);
    end
end


for j = 1:max(clusters)
    ind = clusters == j;
    %     x_final = [x_final; mean(xf(ind))];
    %     y_final = [y_final; mean(yf(ind))];
    %     z_final = [z_final; mean(zf(ind))];
    
    for i = 1:numel(fn)
        if strcmp(fn{i},'fits') ||strcmp(fn{i},'crlbs')
            cdata_out.(fn{i}) = [cdata_out.(fn{i}); mean(cdata_in.(fn{i})(ind,:))];
        else
            cdata_out.(fn{i}) = [cdata_out.(fn{i}), mean(cdata_in.(fn{i})(ind))];
        end
    end
end
cdata_out.xf_loc_unc = std(x_tot);
cdata_out.yf_loc_unc = std(y_tot);
cdata_out.zf_loc_unc = std(z_tot);
end