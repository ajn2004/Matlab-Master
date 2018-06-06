%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Gas Setup
%
% This is a script to get started using neural gases in matlab
% 5/9/16 AJN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% user defined variables
% gas params
nodes = 2; % number of nodes
els0 = 0.1; % starting learning rate
ens0 = 0.001; % starting neighborhood learning rate
agem = 5;  % max age of edges
n0 = 20;
max_dist = 0.5;
max_edge = 20;
% its = 1000000;
% sigma = 80;
% il = 30000;
% in = 300;
lambd = 100; % fixed add node iteration cycle
max_node = 50;
beta = 0.5;
alpha = 0.5;

% 'data' params
points = 10000;
clusts = 3;

%% Data Generation
xf = [];
yf = [];
for i = 1:clusts
    xf = [xf; 10*rand(round(points/clusts),1) + 10*i];
    yf = [yf; (2*rand(round(points/clusts),1)+ 5*i)];
    
end

% xf = [xf; (max(xf)-min(xf))*rand(round(points/3),1) + min(xf)];
% yf = [yf; (max(yf)-min(yf))*rand(round(points/3),1) + min(yf)];

%% Initialize the net
nodex = (max(xf) - min(xf))*rand(nodes,1)+ min(xf);
nodey = (max(yf) - min(yf))*rand(nodes,1)+ min(yf);
edges = [1,2,0]; % edges will be represented by vertex1, vertex2, and age
errors = zeros(nodes,1);
% represented data
plot(xf,yf,'.b')
hold on
plot(nodex, nodey,'.r')
hold off
%% begin clustering

ens = ens0;
els = els0;
i = 0;
while true
    
    inde = randperm(numel(xf));
    xt =  xf(inde(1)); % pick out random vector
    yt = yf(inde(1));
    %% Find minimum distance to random vector
    dists = ((nodex - xt).^2 + (nodey - yt).^2).^0.5;
    inds = find(dists == min(dists)); % index of minimal node s
    indt = find(dists == min(dists(dists>min(dists))));
    indt = indt(1);
    % updating vectors
    errors(inds) = errors(inds) + dists(inds)^2;
    
    dust = ((xf-xt).^2 + (yf - yt).^2).^0.5;
    n = numel(dust(dust < max_dist));
    nodex(inds) = nodex(inds) + els*exp(-n0/n)*(xt - nodex(inds)); % update s position
    nodey(inds) = nodey(inds) + els*exp(-2*n0/n)*(yt - nodey(inds));
    
    % find topological neighbors to s
    edgef = edges(:,1:2);
    [row, col] = find(edgef == inds);
    flag = 0;
    for k = 1:numel(row) % update neighborhood
        for j = 1:2 %find the non-minimal node for updating
            if edges(row(k),j) ~= inds  % if the index is not the one we've looked at already
                indk = edges(row(k),j); % grab index for subsequent use
                nodex(indk) = nodex(indk) + ens*exp(-n0/n)*(xt - nodex(indk)); % update s neighbor position
                nodey(indk) = nodey(indk) + ens*exp(-n0/n)*(yt - nodey(indk));
                if edges(row(k),j) == indt % if we run across second nearest vector t reset age to 0
                    edges(row(k),3) = -1;
                    flag =1;
                end
                edges(row(k),3) = edges(row(k),3) + 1; %increment edge age by 1
            end
        end
    end
    
    if flag == 0 % then s and t were not edges and must be made edges
        edges = [edges;[inds, indt, 0]];
    end
    clear row
    % remove any edges older than age max
    [row] = find(edges(:,3) > agem);
    edges(row,:) = [];
    keepnodes = ones(numel(nodex),1);
    indecies =[];
    edgef = edges(:,1:2);
    for l = 1:numel(keepnodes) % go through each index to see if it has edges, if not remove
        clear row col
        [row, col] = find(edgef == l);
        if isempty(row)
            keepnodes(l) = logical(0);
            indecies(l) = l;
        end
    end
    
    
    nodex(logical(1 - keepnodes)) = [];
    nodey(logical(1 - keepnodes)) = [];
    errors(logical(1 - keepnodes)) = [];
    
    edgemod = edges.*0;
    if ~isempty(indecies)
        for ml = 1:numel(indecies) % go through and count up how many indecies should be dropped for each index
            if indecies(ml) > 0
                for lm = 1:numel(edges(:,1))
                    if edges(lm,1) > indecies(ml)
                        edgemod(lm,1) = edgemod(lm,1) + 1;
                    end
                    if edges(lm,2) > indecies(ml)
                        edgemod(lm,2) = edgemod(lm,2) +1;
                    end
                end
            end
        end
    end
    edges =edges - edgemod;
    if i/lambd == round(i/lambd) && numel(nodex) < max_node % add a node
        inde = find(errors == max(errors));
        % go through nodes to find neighbors connected
        edgef = edges(:,1:2);
        [rows, cols] = find(edgef == inde);
        maxerr2 = -1;
        for k = 1:numel(rows) % search neighborhood
            for j = 1:2 %find the non-minimal node for updating
                if edges(rows(k),j) ~= inde  % if the index is not the one we've looked at already
                    if errors(edges(rows(k),j)) > maxerr2
                        maxerr2 = errors(edges(rows(k),j)); %always pick out maximum error
                        inde2 = edges(rows(k),j);
                        indout = rows(k);
                    end
                end
            end
        end % at this point maxerr2 should contain the second largest error connected to the vector with the largest error
        %         inde2 = find(errors == maxerr2);
        
        % build information of the next point
        nodex = [nodex; (nodex(inde) + nodex(inde2))/2];
        nodey = [nodey; (nodey(inde) + nodey(inde2))/2];
        edges(indout,:) =[];
        edges = [edges;[numel(nodex) , inde, 0]; [numel(nodex), inde2, 0]];
        errors = [errors; errors(inde)];
        errors(inde) = alpha*errors(inde);
        errors(inde2) = alpha*errors(inde2);

    end
    
    els = els0; % * exp(-ic/il);
    ens = ens0; % * exp(-ic/in);
    plot(xf,yf,'.k')
    hold on
    plot(nodex, nodey,'.r')
    for h = 1:numel(edges(:,1))
        plot(nodex(edges(h,1:2)),nodey(edges(h,1:2)),'b');
        
    end
    hold off
    xlabel(['Iteration ', num2str(i),', number of nodes ', num2str(numel(nodex))]);
    drawnow
    %     waitforbuttonpress
    i = i+1;
    % update total error
    errors = errors - beta*errors;
    M(i) = getframe(gcf);
end