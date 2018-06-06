function [clus_id] = func_you_clusters(w, xf_all,yf_all, perc_nn, min_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func_you_clusters
%
% This is a program to make clusters based off of overlapping localizations
%
%
% 8/22/16 AJN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% perc_nn = 0.7;  % what precentage of the average nearest neighbor distance of nodes to use for circle radius in determining member localizations
% min_num = 1; % minimum number of overlapping localizations to be considered apart of the same cluster
% [fname, fpath] = uigetfile('*.mat');
% load([fpath,fname]);

% determine nearest neighbor distribution
for i = 1:numel(w(:,1))
    d2s = ((w(:,1) - w(i,1)).^2 + (w(:,2) - w(i,2)).^2);
    min_NN(i) = min(d2s(d2s>0)).^0.5;
end

dmax = mean(min_NN)*perc_nn;
% clear min_NN

% at this point we have a characteristic distance
clus_id = zeros(numel(xf_all),1);
% now to build clusters based on localization overlaps
for i = 1 : numel(w(:,1))
    
    ds = ((xf_all - w(i,1)).^2 + (yf_all - w(i,2)).^2).^0.5; % determine distances from each point
    
    % find all points within radius dmax
    ind = find(ds <= dmax);
    % at this point we have all the localizations associated with node i
    clus_id_temp = clus_id(ind);
    
    % cluster assignment
    if sum(clus_id_temp) == 0 % if none of the localizations in temp cluster are associated with a cluster assign to cluster max_cluster +1
        clus_id(ind) = max(clus_id) + 1;
        
    else %if the sum of the clus_id_temp is greater than 0 then there begin overlap protocol
        % there are 2 cases of overlap, 1 where there is enough
        % localizations for overalp, and 1 where there is isn't enough
        %         scatter(xf_all,yf_all,10,clus_id)
        
        % first find all clusters identifications
        old_clus_inds = find(clus_id_temp > 0); %gets all clusters above 0
        list_o_clus = unique(clus_id_temp(old_clus_inds)); % this should be a unique list of previously found cluster ids
        ord_list = sort(list_o_clus);  % this is an ordered list
        this_clus = 0; % set the value of this cluster equal to 0 initiially it will serve as a flag for overlapping clusters and a running identification
        change_inds = [];  %this will become the array of clusters that need to be updated because of this node
        
        for j = 1:numel(ord_list) % loop through found cluster ids
            
            cluster = ord_list(j);
            ind_clus = find(clus_id_temp == cluster); %find the number of localizations with relevant identification
            
            if numel(ind_clus) >= min_num % if the number of overlap containing the relevant id, then check to see if another cluster has overlapped
                
                if this_clus == 0 % if this_clus is 0 then no other clusters have been found to overlap with our point i
                    this_clus = cluster;  %update the cluster value to be the first one found (due to the sort command this should ensure the lowest clustered number is chosen
                    change_inds = this_clus;
                    
                else % if this statement is executed it is recognized that another cluster has been found to overlap, and we must adjust that by recording which cluster markings should be changed
                    %given the way we are handling the indexes all we have
                    %to do right now is add the 'new' overlapping cluster
                    %to the list of change_inds
                    change_inds = [change_inds; cluster];
                end
            end
        end
        
        if this_clus ~= 0 % if we found overlapping clusters then we have to rename all clusters who were overlapping
            for j = 1:numel(change_inds)
                ind_chng = change_inds(j);
                indexes = find(clus_id == ind_chng);
                clus_id(indexes) = this_clus;
                clus_id(ind) = this_clus;
            end
            
        else  % if no overlapping clusters were found then set cluster ID to 1 more than previous max
            clus_id(ind) = max(clus_id) + 1;
        end
    end
end

