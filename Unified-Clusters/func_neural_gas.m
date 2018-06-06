function [w, C, T] = func_neural_gas(data, nodes, ei, ef, li, lf, lifeti, lifetf, tmax, graph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func_neural_gas
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




%% Following Neural Gas Paper
N = numel(data(:,1));
dim = numel(data(1,:));

%% Step 0: 
%Assign Initial values to the weights wi belonging to R^n and set
% all Cij to 0
fig = figure('units','normalized','outerposition',[0 0 1 1],'Visible','Off');
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
count = 1;
for j = 1:dim % loop over the dimensions of the data
    w(:,j) = (max(data(:,j)) - min(data(:,j)))*rand(nodes,1)+ min(data(:,j));
end

t = 0; % start of time
dex = 1:9;
for i = 1:4
dex = [dex,10^(i):0.1*10^(i):9*10^i];
end
dex(dex>tmax) = [];
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
   if C(i0,i1) == 1
       T(i0,i1) = 0;
       T(i1,i0) = 0;
   end
   if C(i0,i1) == 0
       C(i0,i1) = 1;
       C(i1,i0) = 1;
       T(i0,i1) = 0;
       T(i1,i0) = 0;
   end
  
   %% Step 5
   % Increase the age of all connections of i0 by incrementing the i0th row
   % of the T matrix by 1
   clear cols rows
   [cols] = find(C(i0,:) > 0); % this gives the connections of the nodes
   [rows] = find(C(:,i0) > 0); % this gives the connections of the nodes
   
   if ~isempty(cols) && sum(cols > numel(C(1,:))) < 1
       for ii = 1:numel(cols)
           if cols(ii) <= numel(C(1,:)) && cols(ii) > 0  % Ensure the row / col variable stays within range
               T(i0,cols(ii)) = T(i0,cols(ii)) + 1; % increment each connection age by 1
               T(rows(ii),i0) = T(rows(ii),i0) + 1; % increment each connection age by 1
           end
       end
   end
   clear cols rows
   %% Step 6
   % Remove all connections of i0 which exceeded their lifetime by setting
   % C(i0,j) = 0 for all j with C(i0,j) = 1 and T(i0,j) > lifet.
   [cols] = find(C(i0,:) == 1 & T(i0,:) > lifet);
   [rows] = find(C(:,i0) == 1 & T(:,i0) > lifet);
   C(i0,cols) = 0;
   C(rows,i0) = 0;
   
   if sum(dex == t) == 1
   % This will be used only for creating the figure for munc13
   s = plot(data(:,1),data(:,2),'.k');
   hold on
   plot(w(:,1),w(:,2),'.r')
   title(['Neural Gas after ', num2str(t),' iterations']);
   [rows, cols] = find(C > 0); % here we plot connection grid
   for i = 1:numel(rows)
       plot([w(rows(i),1),w(cols(i),1)],[w(rows(i),2),w(cols(i),2)],'b','LineWidth',1.5);
   end
   plot(w(:,1),w(:,2),'.r','markersize',15); % plot positions of nodes
   xlabel(['Using ', num2str(nodes),' nodes']);
   hold off
   legend('Data','Nodes','Connections');
%    M(count) = getframe(fig);
   count = count +1;
   end
    
   t = t+1;
   
   if t > tmax
       break
   end
end
% movie2gif(M,'makeadamovie.gif','DelayTime',0.02);
% if strcmp(graph,'y') || strcmp(graph,'Y')
%     fig = figure('units','normalized','outerposition',[0 0 1 1]);
%     if dim == 2
%        plot(data(:,1),data(:,2),'.k'); % plot data
%        hold on
%        plot(w(:,1),w(:,2),'.r'); % plot positions of nodes
%        title(['Final Result of Neural Gas after ', num2str(t),' iterations']);
% %        t
% %        %    waitforbuttonpress
%        [rows, cols] = find(C > 0); % here we plot connection grid
%        for i = 1:numel(rows)
%            plot([w(rows(i),1),w(cols(i),1)],[w(rows(i),2),w(cols(i),2)],'b','LineWidth',1.5);
%        end
%        plot(w(:,1),w(:,2),'.r','markersize',15); % plot positions of nodes
%        xlabel(['Using ', num2str(nodes),' nodes']);
%        hold off
%        legend('Data','Nodes','Connections')
%        drawnow
%     elseif dim >=3
%         scatter3(data(:,1),data(:,2),data(:,3),'k');
%         hold on
%         
%     end
% end



end
