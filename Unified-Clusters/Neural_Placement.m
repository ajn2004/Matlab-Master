%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Placing
%
% This is script that will train a neural capable of placing data vectors
% into the clustered regions mapped out in GeNeural Cluster
%
%
%
%
% AJN 5/12/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%% User Controlled Variables
lambda = 5;
trainper = 0.8;
epsint = 0.12;
hidden_layer_size = 30;

[fname, fpath] = uigetfile('*clustered*.mat');

cd(fpath);
load(fname);

%% Data Organization
train_data = [w(:,1), w(:,2)];
for i = 0.1:0.1:8
    train_data = [train_data,sin(2*pi*w(:,1)/(i)), cos(2*pi*w(:,1)/(i)),sin(2*pi*w(:,2)/(i)), cos(2*pi*w(:,2)/(i))];
end
nodes = numel(train_data(:,1));

switchdex = randperm(nodes);

trs = round(nodes*trainper);

train_data0 = train_data(switchdex(1:trs),:);
train_data1 = train_data(switchdex(trs + 1:end),:);
y0 = clus_id(switchdex(1:trs));
y1 = clus_id(switchdex(trs +1 :end));
clus_id0 = clus_id(switchdex(1:trs));
clus_id1= clus_id(switchdex(trs+1:end));
input_layer_size = numel(train_data(1,:));

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

output_layer_size = numel(index)+1;


% build random thetas
Theta1 = rand(hidden_layer_size, input_layer_size+1)*2*epsint - epsint;
Theta2 = rand(output_layer_size, hidden_layer_size+1)*2*epsint - epsint;

thetas0 = [Theta1(:);Theta2(:);];
options = optimset('MaxIter', 1000);
[thetasf, cost ] = fmincg(@(p)neural_multi_class_cost(p,input_layer_size,hidden_layer_size,output_layer_size, train_data0, clus_id0, index, lambda), thetas0, options);

Theta1 = reshape(thetasf(1:(input_layer_size + 1)*hidden_layer_size), hidden_layer_size, input_layer_size + 1);
Theta2 = reshape(thetasf((input_layer_size + 1)*hidden_layer_size +1 : end), output_layer_size, hidden_layer_size + 1);

% test on training data and testing data
[train_clus, trh2] = neural_predict(Theta1, Theta2, train_data0);
[test_clus, teh2] =  neural_predict(Theta1, Theta2, train_data1);


% All non clusters will be represented as 0 cluster and must be updated to
% compare

index = [index,0];
indout = find(counts <= ask1);
indo = [];
for i = 1:numel(indout)
        indo = [indo;find(clus_id == indout(i))];
end

clus_id_lim = clus_id;
clus_id_lim(indo) = 0;
y0 = clus_id_lim(switchdex(1:trs));
y1 = clus_id_lim(switchdex(trs +1 :end));
disp(['Accuracy on training set is ', num2str(100*mean(double(y0 == index(train_clus).'))), '%']);
disp(['Accuracy on testing set is ', num2str(100*mean(double(y1 == index(test_clus).'))), '%']);
