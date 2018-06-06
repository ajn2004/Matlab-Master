% function [Theta1n, Theta2n] = teach_this_net(Theta1, Theta2, train_data, index ,clus_id, trainper, lambda)
    % Prepare the data
     X = [ones(numel(train_data(:,1)),1),train_data];
     Y = zeros(numel(train_data(:,1)),numel(Theta2(:,1)));
     


% end