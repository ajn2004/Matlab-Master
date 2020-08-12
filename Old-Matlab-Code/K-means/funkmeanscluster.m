function [clus_coords, allclosest, M] = funkmeanscluster(noc, data0, showdata, viewvec)
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

%% Argument management
if nargin < 2
    error('Too few arguments, make sure you input (number of clusters, data) at a minimum');
elseif nargin >4
    error('Too many arguments, make sure you input (number of clusters, data, showdata, view angle) at maximum')
end

if nargin == 2
    showdata = 'n';
    viewvec = [0 0];
end

if nargin == 3
    viewvec = [23,21];
end
%% Get information about the data set
ac= noc;
numvar = numel(data0(1,:));
m = numel(data0(:,1));
clus_check = zeros(noc,numvar);
check = 0;
count = 1;
% if numvar < 2
%     showdata = 'n';
% end
if m > 10000
    mm = 10000;
else
    mm = m;
end
fprintf('Crunching');
p = randperm(m);
%% generate initial guess for cluster locations based on data
if strcmp(showdata,'y') || strcmp(showdata, 'Y')
    figure
end
while true
    
    
    data = data0(p(1:mm),:);
    q = randperm(mm);
    for i = 1:noc
        for j = 1:numvar
            clus_coords(i,j) = data0(q(i),j);
        end
    end
    
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
        for i = 1:mm
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
            if strcmp(showdata,'y') || strcmp(showdata,'Y')
                if numvar > 3
                    showvar = 3;
                elseif numvar == 2
                    showvar = 2;
                else
                    showvar = 1;
                end
                switch showvar
                    case 3
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
                    case 2
                        if i == 1
                            plot(data(index,1), data(index,2), '.g');
                            hold on
                        elseif i == 2
                            plot(data(index,1), data(index,2), '.b');
                        elseif i == 3
                            plot(data(index,1), data(index,2), '.c');
                        elseif i == 4
                            plot(data(index,1), data(index,2), '.m');
                        end
                    case 1
                        if i == 1
                            plot(data(index,1), '.g');
                            hold on
                        elseif i == 2
                            plot(data(index,1), '.b');
                        elseif i == 3
                            plot(data(index,1), '.c');
                        end
                    otherwise
                end
            end
            
        end
        if strcmp(showdata,'y') || strcmp(showdata,'Y')
            switch showvar
                case 3
                    scatter3(last_guess(:,1),last_guess(:,2), last_guess(:,3), '.y');
                    scatter3(clus_coords(:,1), clus_coords(:,2), clus_coords(:,3), '.r');
                    view(viewvec);
                case 2
                    plot(last_guess(:,1),last_guess(:,2), '.y');
                    plot(clus_coords(:,1), clus_coords(:,2), '.r');
                otherwise
            end
            hold off
            axis image
            
            drawnow
            M(count) = getframe(gcf);
        end
        if isequal(last_guess,clus_coords)
            break
        end
        count = count+1;
        fprintf('.');
        if count/60 == round(count/60)
            fprintf('\n');
        end
    end
    
    %% reverify answers in cluster by checking against previous results
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
fprintf('\n')
%% build index for final output
clear distmat
for i = 1:noc
    for j = 1:numvar
        distmat(:,j,i) = (data0(:,j) - clus_coords(i,j)).^2;
    end
end

clear distance
for j = 1:noc
    %% PICK UP HERE
    distance(:,j) = sum(distmat(:,:,j),2).^0.5;
end
[c, allclosest] = min(distance,[],2);
disp('Finished');
% for i = 1:noc
%     disp(['% of elements in cluster ', num2str(i),' is ', num2str(mean(allclosest == i))]);
% end
end