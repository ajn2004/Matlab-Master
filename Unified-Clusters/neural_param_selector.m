%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural_gas param selection
%
% This is a matlab script that attempts to exactly recreate the algorithm
% described in Martinetz and Schulten 1991
%
% This version of the script has been turned into a function 5/11/16 AJN
%
% Data should be organized in an data(dataindex,dimension index) matrix format
%
% AJN 5/10/16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nodes = 3000; % The number of nodes making up your gas  reasonable: 800
epsilon_i = 0.3; % step size in weight vector learning eq 1. should be between 0 and 1
epsilon_f = 0.01;  % reasonable values are e_i = 0.3 and e_f = 0.01
                    % This variable controls the extent of the step size,
                    % too large inital variable will result into
                    % unreasonable clumping of many of the nodes

lambdai = 1000;  % decay rate vector for exponential decay distance updating
lambdaf = 0.01;  % This variable affects how much the nodes will expand into the less dense regions
                % If the initial value is too low relative to your number
                % of nodes, you will see more nodes left outside the region
                % of interest and not connected at all
                % Too high final values will not allow for proper filling of
                % clustered strutcure. Resonable values lambdai = 0.5*nodes,
                % lambdaf = 0.5

lifeti = 20;  % maximum lifetime of any connection without being refreshed
lifetf = 3;   % This variable controls how connected the web is in the end, if these values are too high
             % You will notice each node having too many connections, if
             % too low you will not have connected structure
               % reasonable values are lifeti = 20 lifetf = 80

tmax = 40000;   % this is the final iteration  reasonable: 40000

names = {'xf_all', 'yf_all','q'}; % names of variables

% user Independent
[fname, fpath] = uigetfile('*tol.mat');
addpath(pwd);
cd(fpath)
    close all
    for i = 1:numel(names)
        load(fname,names{i});
    end
    data = [xf_all*q,yf_all*q];
    

%% Following Neural Gas Paper
N = numel(data(:,1));
dim = numel(data(1,:));

%% Step 0: 
%Assign Initial values to the weights wi belonging to R^n and set
% all Cij to 0

% Notes to self
% In this case we are doing a 2 dimensional case of [x position, and
% y position] for a weighting vector where each 'node' is a vector defined
% by us.
% we will also be only using part of the C matrix, which will continue to
% build

k = 0:nodes-1;    % initializing k vector which determines the number of nodes an index i has have a distance less than it
w = zeros(nodes,2);% weight vector initialization
C = zeros(nodes); % C vector is set to 0s here we examine connections between [ith vector, jth vector]
T = C;            % T vector tracks the age of each connection
distances = k;
% randomly assign weight vectors 

for j = 1:dim % loop over the dimensions of the data
    w(:,j) = (max(data(:,j)) - min(data(:,j)))*rand(nodes,1)+ min(data(:,j));
end

t = 0; % start of time

%% Begin Neural Gas Iteration Algorithm
while true 
    %% prestep 1 updating variables
    epsilon = (epsilon_i*(epsilon_f/epsilon_i)^(t/tmax));
    lambda = (lambdai*(lambdaf/lambdai)^(t/tmax));
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
   
   t = t+1;
   if t > tmax
       break
   end
end


    figure('units','normalized','outerposition',[0 0 1 1]);
    if dim == 2
       plot(data(:,1),data(:,2),'.k'); % plot data
       hold on
       plot(w(:,1),w(:,2),'.r'); % plot positions of nodes
       title(['Final Result of Neural Gas after ', num2str(t),' iterations']);
%        t
%        %    waitforbuttonpress
       [rows, cols] = find(C > 0); % here we plot connection grid
       for i = 1:numel(rows)
           plot([w(rows(i),1),w(cols(i),1)],[w(rows(i),2),w(cols(i),2)],'b','LineWidth',1.5);
       end
       plot(w(:,1),w(:,2),'.r','markersize',15); % plot positions of nodes
       xlabel(['Using ', num2str(nodes),' nodes']);
       hold off
       legend('Data','Nodes','Connections')
       drawnow
    elseif dim >=3
        scatter3(data(:,1),data(:,2),data(:,3),'k');
        hold on
        
    end
    end
