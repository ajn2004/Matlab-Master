function cost = nn_drift_cost(sub_locs_1,sub_locs_2)
% This will determine the cost of two localization sub_groups using the
% probability based nearest neighbor distances
% bins = 100;
% a = 100;
% combo_sub = [sub_locs_1;sub_locs_2];
% % standard_x = mod(a,a*(combo_sub(:,1) - min(combo_sub(:,1)))/(max(combo_sub(:,1))- min(combo_sub(:,1))));
% % standard_y = mod(a,a*(combo_sub(:,2) - min(combo_sub(:,2)))/(max(combo_sub(:,2))- min(combo_sub(:,2))));
% standard_x = combo_sub(:,1);
% standard_y = combo_sub(:,2);
% % standard_z = mod(a,a*(combo_sub(:,1) - min(combo_sub(:,1)))/(max(combo_sub(:,1))- min(combo_sub(:,1))));
% [f_x] = histcounts(standard_x,bins);
% [f_y] = histcounts(standard_y,bins);
% % histogram(standard_y,bins)
% 
% prob_x = f_x./sum(f_x);
% prob_y = f_y./sum(f_y);
% ind_x = prob_x >0;
% ind_y = prob_y >0;
% cost = -sum(prob_x(ind_x).*log(prob_x(ind_x))) - sum(prob_y(ind_y).*log(prob_y(ind_y)));



bins = 10;
[~,nn_dist] = knnsearch(sub_locs_2,sub_locs_1,'K',4);
nn_dist = mean(nn_dist,2);
% nn_dist(end+1) = 200;
% create bin uniform bins to determine the weighting of the nn distances
bin_width = (max(nn_dist) - min(nn_dist))/(bins - 1);
frequency(1:bins,1) = 0;
which_bin = nn_dist*0;
est_bins = (nn_dist - min(nn_dist))/bin_width;
for i = 1:numel(nn_dist) % populate the histogram and keep track of where each entry bins to
%     if round(est_bins(i))+1 < bins
        frequency(round(est_bins(i))+1) = frequency(round(est_bins(i))+1)+1;
        which_bin(i) = round(est_bins(i))+1;
%     else
%         frequency(end) = frequency(end) + 1;
%         which_bin(i) = bins;
%     end
end
    probability = frequency./sum(frequency);
%     cost = sum(probability(which_bin).*nn_dist);
cost = mean(nn_dist);
%     cost = sum(probability.*log(probability));
end

