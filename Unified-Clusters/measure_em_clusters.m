%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure'em clusters
%
%
% This is a script to measure clusters identified through Unified Clusters
% We'll need to know number of localizations, area, perimeter
%
% AJN 9-8-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;


%% Choose and load files
[fname, fpath] = uigetfile('*clust.mat');
load([fpath,fname]);
data = double([xf_all*q, yf_all*q]);
%% Loop over clusters
for i = 1:max(clus_id)
    % grab the index of all clusters with id i
    index = clus_id == i;
    sub_data = [ data(index,1), data(index,2)];
    c(i).sub_data = sub_data;
    % get number of localizations in cluster
    mem_num(i,1) = numel(sub_data(:,1));
    [rotation{i}, area(i,1)] = boundary(sub_data(:,1),sub_data(:,2));
    c(i).area = area(i,1);
    c(i).mem_num = mem_num(i,1);
    % compute the perimeter from the rotation coords
    perim = 0;
    
    for j = 1:numel(rotation{i})-1
        c(i).orbit(j,:) = [sub_data(rotation{i}(j),1),sub_data(rotation{i}(j),2)];
        dist = ((sub_data(rotation{i}(j),1) - sub_data(rotation{i}(j+1),1)).^2 + (sub_data(rotation{i}(j),2) - sub_data(rotation{i}(j+1),2)).^2).^0.5;
        perim = perim + dist;
    end
    if ~isempty(j)
        c(i).orbit(j+1,:) = [sub_data(rotation{i}(1),1),sub_data(rotation{i}(1),2)];
    end
    perimeter(i,1) = perim;
    c(i).perimeter = perimeter(i,1);
end

for j = 1:1000
    index = area > j/1000;
    num(j) = sum(index);
end

for j = 1:1000
    index = mem_num > j;
    num_mem(j) = sum(index);
end

for j = 1:10000
    index = perimeter > j/1000;
    num_per(j) = sum(index);
end

subplot(1,3,1);plot((1:1000),num_mem)
xlabel('Minimum number of localizations')
ylabel('Number of clusters remaining');
title('Minimum Number of Localizations');
subplot(1,3,2);plot((1:1000)./1000,num)
xlabel('Area in um^2')
ylabel('Number of clusters remaining');
title('Minimum Area');
subplot(1,3,3);plot((1:10000)./1000,num_per)
xlabel('Perimeter in um')
ylabel('Number of clusters remaining');
title('Minimum Perimeter');

clear dist perim rotation sub_data index
save([fname(1:end-4),'_meas.mat'],'c')
clearvars -except c