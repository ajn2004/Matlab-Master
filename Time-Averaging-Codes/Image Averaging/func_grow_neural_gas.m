function [cents] = func_grow_neural_gas(im1)
%% Data Preprocessing
midval = mean(im1(:));
count = 1;
for i = 1:numel(im1(:,1))
    for j = 1:numel(im1(1,:))
        if im1(i,j) > midval
            pos(count,:) = [j,i];
            count = count+1;
        end
    end
end

% At this point pos should represent the positions of pixels that are
% bigger than the average value of im1, this should
%% Begin Neural Gas Iteration Algorithm
while true 
    %% prestep 1 updating variables
    epsilon = (ei*(ef/ei)^(t/tmax));
    lambda = (li*(lf/li)^(t/tmax));
    lifet = (lifeti*(lifetf/lifeti)^(t/tmax));
    %% Step 1
    % Select an input vector v of the input manifuld M.
    
    % Notes to self
    % Select a random data vector
    index = randi(numel(data(:,1)));
    data0 = data(index,:);
    
    %% Step 2
    % For each unit i determine the number of ki neural nets with distances
    % smaller than their
    % build a vector of distances where each row corresponds to the
    % distance from the ith node
    
   
    for i = 1:nodes % determine the euclidean distance for each node to the data point
        ds2 = 0;
        for j = 1:dim
            ds2 = ds2 + (data0(j) - w(i,j))^2;
        end
        distances(i) = ds2^0.5;
    end

    % get ranking of distances
    [rank_dist, backrank, ranking] = unique(distances);
    k = ranking - 1; % the ranking vector orders from smallest distance to least, with the smallest being 1. Subtracting each value by 1 tells how many nodes are closer than node i.
    
   %% Step 3
   % Perform an adaptation step for the weights
   for i = 1:nodes
        for j = 1:dim
            w(i,j) = w(i,j) + epsilon*exp(-k(i)./lambda).*(data0(j) - w(i,j));
        end
   end
   %% Step 4
   % If C(i0,i1) =0 set to 1, and T(io,i1) = 0 If C(i0,i1) = 1 set T(i0,i1) = 0
   i0 = find(k == 0); % grab minimum index
   i1 = find(k == 1); % grab second minmum index 
   if C(i0,i1) == 1;
       T(i0,i1) = 0;
       T(i1,i0) = 0;
   end
   if C(i0,i1) == 0;
       C(i0,i1) = 1;
       C(i1,i0) = 1;
       T(i0,i1) = 0;
       T(i1,i0) = 0;
   end
  
   %% Step 5
   % Increase the age of all connections of i0 by incrementing the i0th row
   % of the T matrix by 1
   [cols] = find(C(i0,:) > 0); % this gives the connections of the nodes
   [rows] = find(C(:,i0) > 0); % this gives the connections of the nodes
   
   T(i0,cols) = T(i0,cols) + 1; % increment each connection age by 1
   T(rows,i0) = T(rows,i0) + 1; % increment each connection age by 1
   clear cols
   %% Step 6
   % Remove all connections of i0 which exceeded their lifetime by setting
   % C(i0,j) = 0 for all j with C(i0,j) = 1 and T(i0,j) > lifet.
   [cols] = find(C(i0,:) > 0 & T(i0,:) > lifet);
   [rows] = find(C(:,i0) > 0 & T(:,i0) > lifet);
   C(i0,cols) = 0;
   C(rows,i0) = 0;
   
   
   t = t+1;
   if t > tmax
       break
   end
end

end