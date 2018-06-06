%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVM Placing
%
% This is script that will train a support vector machine capable of placing data vectors
% into the clustered regions mapped out in GeNeural Cluster
%
%
%
%
% AJN 5/16/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%% User Controlled Variables

trainper = 0.8;


[fname, fpath] = uigetfile('*clustered*.mat');

cd(fpath);
load(fname);

%% Data Organization
train_data = [w(:,1), w(:,2)];
nodes = numel(train_data(:,1));
switchdex = randperm(nodes);
trs = round(nodes*trainper);
train_data0 = train_data(switchdex(1:trs),:);
train_data1 = train_data(switchdex(trs + 1:end),:);
y0 = clus_id(switchdex(1:trs));
y1 = clus_id(switchdex(trs +1 :end));
clus_id0 = clus_id(switchdex(1:trs));
clus_id1= clus_id(switchdex(trs+1:end));

%% Find the number of molecules in each cluster and throw out those that have too few members
for i = 1:max(clus_id)
    counts(i) = sum(clus_id == i);
end


% in a sense the training vectors are already built and even segmented All
% we have to do is use them to train

%% Cluster Cleanup

% Show changes in average counts is related to cutoff distance
for i = 1 : 80
    y(i) = mean(counts(counts>i));
end


while true
    ask1 = input('What cutoff size would you like to try? ');
    index = find(counts > ask1);
    subplot(1,2,1);plot(1:80,y);
    title('Mean number of nodes / cluster vs. Cutoff distance');
    xlabel('Cutoff Distance nm');
    ylabel('Average cluster size');
    hold on
    subplot(1,2,1); plot([ask1,ask1],[0, max(y)],'r');
    disp(['There are ', num2str(numel(index)),' clusters at that cutoff with average weight',mean(counts(counts > ask1)), ' molecules']);
    
    hold off
    
    ind = [];
    for i = 1:numel(index)
        ind = [ind;find(clus_id == index(i))];
    end
    subplot(1,2,2); plot(xf_all*q,yf_all*q,'.k');
    subplot(1,2,2); scatter(w(ind,1),w(ind,2),15,clus_id(ind),'Filled');
    ask2 = input('Are you happy with this distance? ', 's');
    if strcmp(ask2,'y') || strcmp(ask2,'Y')
        break; % at this point we have an index of the clusters already arranged for us
    end
end

% At this point we have ind which gives us every vector, and the cluster it
% is associated with train_data(ind,:) corresponds to cluster
% clus_id(ind)

clus_num = numel(index)+1;
clear ind
%% 1 v All SVM Network
% The principle is that we will use matlab software to train a support
% vector network to separate clus_num areas in the data
m = numel(train_data0(:,1)); % number of examples
n = clus_num;
Y = zeros(m,n);
Y(:,end) = 1; % last component is a "trash" component

for i = 1:n-1 % convert vectors into 0 and 1s and train column's SVMmodel
    ind = y0 == index(i);
    Y(ind,i) = 1;
    Y(ind,end) = 0; 
    eval(['[SVMModel_',num2str(i),'] = fitcsvm(train_data0,Y(:,i),','''','KernelFunction', '''',',','''','rbf','''',');'])
end

% Determination of effectiveness of the SVM network on remaining training set
m = numel(train_data1(:,1)); % number of examples
Y = zeros(m,n);
Y(:,end) = 1; % last component is a "trash" component
labels = [];
for i = 1:n-1 % convert vectors into 0 and 1s and train column's SVMmodel
    ind = y1 == index(i);
    Y(ind,i) = 1;
    Y(ind,end) = 0; 
    eval(['[label, score] = predict(SVMModel_',num2str(i),',train_data1);'])
    labels = [labels,label];
end

% translate labels to cluster list
for i = 1:m
    [val, ind] = max(labels(i,:));
    [valy, indy] = max(Y(i,:));
    if val > 0.2
        p(i) = ind;
        py(i) = indy;
    else
        p(i) = n;
        py(i) = indy;
    end
end

disp(['When measuring the accuracy of placing the SVM network, the percentage is ', num2str(100*mean(p == py)),'%'])
% compare cluster list to actual clusters



% All non clusters will be represented as 0 cluster and must be updated to
% compare

index = [index,0];
indout = find(counts <= ask1);
indo = [];
for i = 1:numel(indout)
    indo = [indo;find(clus_id == indout(i))];
end

data = [xf_all*q, yf_all*q];

% Determination of effectiveness of the SVM network on remaining training set
m = numel(data(:,1)); % number of examples
Y = zeros(m,n);
labels = [];
for i = 1:n-1 % convert vectors into 0 and 1s and train column's SVMmodel
    eval(['[label, score] = predict(SVMModel_',num2str(i),',data);'])
    labels = [labels,label];
end

% translate labels to cluster list
for i = 1:m
    [val, ind] = max(labels(i,:));
    if val > 0.2
        clusters(i,1) = ind;
    else
        clusters(i,1) = n;
    end
end

figure
scatter(xf_all,yf_all,5,clusters,'Filled');