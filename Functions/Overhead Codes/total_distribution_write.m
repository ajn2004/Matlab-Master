function total_distribution = total_distribution_write(total_distribution, red_dist, orange_dist, normals, cluster, clusts)
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
    regions = file_id*0 + max(clusts); % keep track of number of regions
    region_id = file_id*0 + cluster;  % keep track of region
     
    % instead of 'assigning' we have to 'add'
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