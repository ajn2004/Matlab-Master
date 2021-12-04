function update_total_distribution(red_dist, orange_dist, normals, cluster, max_clust, file_id)
% A writing function to keep the composite data object clean as updates to
% algorithms go on. This function call will load the data object
% total_distribution, clear any previous entries for the file/region combo,
% and then append the provided data [red_dist, orange_dist, normals, cluster, max_clust] 
index_file_path = '..\..\Composites\';
index_file_name = 'gpi_vGlut_distribution.mat';
load([index_file_path,index_file_name])
% we want to start with a complete wipe of the previously written 

colors = {'red', 'orange'};
% This will loop over colors and utilize string eval functions
for i  = 1:2
    % find previous entries
    ind = find(total_distribution.(colors{i}).file_id == file_id_number & total_distribution.(colors{i}).region_id == region);
    fnames = fieldnames(total_distribution.(colors{i}));
    for f = 1:numel(fnames)
        if f == 2
            total_distribution.(colors{i}).(fnames{f})(ind,:) = [];
        else
            total_distribution.(colors{i}).(fnames{f})(ind) = [];
        end
    end
    
    % Prepare 'write' process
    str0 = ['file_id = ', colors{i},'_dist(:,1)*0 + ', num2str(file_id_number),';'];
    eval(str0); % build and eval file_id array and define region arrays
    regions = file_id*0 + max_clust; % keep track of number of regions
    region_id = file_id*0 + cluster;  % keep track of region
     
    % This could probably be re_written, but handling strings allowed the
    % iteration through different variables when this code was in
    % live-script stage, so that's what is happening here. This is just a
    % series of appending statements
    str1 = ['total_distribution.',colors{i},'.distances = [total_distribution.',colors{i},'.distances;',colors{i},'_dist(:,1)];'];
    eval(str1);
    str2 = ['total_distribution.',colors{i},'.regions = [total_distribution.',colors{i},'.regions; regions];'];
    eval(str2);    
    str3 = ['total_distribution.',colors{i},'.normals = [total_distribution.',colors{i},'.normals; normals(',colors{i},'_dist(:,2),:)];'];
    eval(str3);    
    str4 = ['total_distribution.',colors{i},'.region_id = [total_distribution.',colors{i},'.region_id; region_id];'];
    eval(str4);
    str5 = ['total_distribution.',colors{i},'.file_id = [total_distribution.',colors{i},'.file_id; file_id];'];
    eval(str5);
    times = file_id*0 + now;
    total_distribution.(colors{i}).times = [total_distribution.(colors{i}).times; times]; % keep track of when this data point was added
end
save([index_file_path,index_file_name])

% 'Return' section
assignin('base','total_distribution',total_distribution);

%% Original live script code
% At this point both red_dist and orange_dist are 'peak normalized'
% We're going to save them
% We need to develop an id system to track which files have been logged and
% not
% To do this we need to agree upon a system
% this distance output is [dist, idx] where idx corresponds to the index
% for the proper center and normal vector associated w/ the data point
% so, our data structure will be
% [data,color,center_vector(1:3),normalvector(1:3), region_id, regions,file_id]

% We are interesting in measuring distribution of distances from cell
% surface, and we also have orientation to surface, so saving distances and
% normal vectors will be the interesting values to begin with. As the
% measurements are meant to bring us into a 'relative space' and making
% reasonable storage decisions are important

% Saving
% red_dist, orange_dist, red_normals, orange_normals, file_info,
% region_info
% 
% Total_distribution Object
% total_distribution.{red, orange}.{distances, normals, regions, file_id, date_written}
% index_file_path = '..\..\Composites\';
% index_file_name = 'gpi_vGlut_distribution.mat';
% load([index_file_path,index_file_name])
% % we want to start with a complete wipe of the previously written 
% 
% colors = {'red', 'orange'};
% % This will loop over colors and utilize string eval functions
% for i  = 1:2
%     % find previous entries
%     ind = find(total_distribution.(colors{i}).file_id == file_id_number & total_distribution.(colors{i}).region_id == region);
%     fnames = fieldnames(total_distribution.(colors{i}));
%     for f = 1:numel(fnames)
%         if f == 2
%             total_distribution.(colors{i}).(fnames{f})(ind,:) = [];
%         else
%             total_distribution.(colors{i}).(fnames{f})(ind) = [];
%         end
%     end
%     
%     % Prepare 'write' process
%     str0 = ['file_id = ', colors{i},'_dist(:,1)*0 + ', num2str(file_id_number),';'];
%     eval(str0); % build and eval file_id array and define region arrays
%     regions = file_id*0 + max(clusts); % keep track of number of regions
%     region_id = file_id*0 + cluster;  % keep track of region
%      
%     % instead of 'assigning' we have to 'add'
%     str1 = ['total_distribution.',colors{i},'.distances = [total_distribution.',colors{i},'.distances;',colors{i},'_dist(:,1)];'];
%     eval(str1);
%     str2 = ['total_distribution.',colors{i},'.regions = [total_distribution.',colors{i},'.regions; regions];'];
%     eval(str2);
%     str3 = ['total_distribution.',colors{i},'.normals = [total_distribution.',colors{i},'.normals; normals(',colors{i},'_dist(:,2),:)];'];
%     eval(str3);    
%     str4 = ['total_distribution.',colors{i},'.region_id = [total_distribution.',colors{i},'.region_id; region_id];'];
%     eval(str4);
%     str5 = ['total_distribution.',colors{i},'.file_id = [total_distribution.',colors{i},'.file_id; file_id];'];
%     eval(str5);
%     times = file_id*0 + now;
%     total_distribution.(colors{i}).times = [total_distribution.(colors{i}).times; times]; % keep track of when this data point was added
% end
% save([index_file_path,index_file_name])