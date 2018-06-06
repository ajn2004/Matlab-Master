% function [coords] = kmeanscluster(noc, data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K means clustering
%
% This is a simplistic algorithm for unsupervised learning
% This will take an arbitrary number of variables and a user specified
% number of intended clusters to independently determine clustering of data
% points
%
% AJN 11/20/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

number_of_clusters = 2; % This is the number of clustering regions you intend to find in the data

noc = number_of_clusters;
clear number_of_clusters

ac= noc;
points = 1000;

%% For testing Purposes
data = [[],[]];
numvar = 3;
build_vec = [];
for i = 1:ac
    for j = 1:numvar
        x = randn(points,1) + 20*rand*i;
        build_vec = [build_vec, x];
    end
    data= vertcat(data, build_vec);
    build_vec = [];
end

m = numel(data(:,1));
clus_check = zeros(noc,numvar);
check = 0;
count = 1;
%% generate initial guess for cluster locations based on data
while true
    p = randperm(m);
    for i = 1:noc
        for j = 1:numvar
            clus_coords(i,j) = data(p(i),j);
        end
    end
    
    

    %
    while true
        last_guess = clus_coords; % store last guess of cluster coords
        % build matrix of squared differences
        for i = 1:noc
            for j = 1:numvar
                distmat(:,j,i) = (data(:,j) - clus_coords(i,j)).^2;
            end
        end
        
        % calculate distances in variable space and return which cluster each
        % data point is closest to
        for i = 1:m
            for j = 1:noc
                distance(i,j) = sum(distmat(i,:,j))^0.5;
            end
            [c, ind] = min(distance(i,:));
            closest(i,1) = ind;
        end
        
        for i = 1:noc
            index = closest == i;
            for j = 1:numvar
                cm(1,j) = sum(data(index,j))/numel(data(index,j));
            end
            clus_coords(i,:) = cm;
            if i == 1
                scatter3(data(index,1), data(index,2),data(index,3), 'g');
                hold on
            elseif i == 2
                scatter3(data(index,1), data(index,2),data(index,3), 'b');
            elseif i == 3
                scatter3(data(index,1), data(index,2),data(index,3), 'c');
            elseif i == 4
                scatter3(data(index,1), data(index,2),data(index,3), 'm');
            end
            
        end
        scatter3(last_guess(:,1),last_guess(:,2), last_guess(:,3), '.y');
        scatter3(clus_coords(:,1), clus_coords(:,2), clus_coords(:,3), '.r');
        
        hold off
        axis image
        view(23,21)
        M(count) = getframe(gcf);
        if isequal(last_guess,clus_coords)
            break
        end
        count = count+1;
    end
    
    %% check the cluster algorithm didn't get stuck
    for i = 1:noc
        flag = 1;
        for j = 1:numel(clus_check(1,1,:))
            if 1-ismember(round(clus_coords(i,:)),round(clus_check(:,:,j)),'rows')
                flag = 0;
            end
        end
        
    end
    if flag == 1
        break
    else
        check = check +1;
        clus_check(:,:,check) = clus_coords;
    end
    
end