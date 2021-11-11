function [dense_data] = add_average_points(data, iters)
%  This is a function to artificially add points to outer hulls
%   Detailed explanation goes here
k = boundary(data,1); % Get the boundary points of the data distribution

% Repeate for iters iterations
for j = 1:iters
    new_points = []; % preallocate a new 
    
    for i = 1:numel(k(:,1))
        new_points = [new_points; mean(sub_data(k(i,:),:))];
    end
    sub_data = [sub_data; new_points];
    k = boundary(sub_data,1);
end
dense_data = sub_data;
end

