function [clust] = DB_scan(data, eps, minp)
% This function applies the DB scan algorithm

C = 0;

clust = (1:numel(data(:,1)))*0;

% scatter(data(:,1),data(:,2),2,[0 0 0],'Filled')
count = 1;

for i = 1:numel(data(:,1))
    if clust(i) == 0 % If we haven't visited the point yet, perform analysis
        N = func_range_scan(data(i,:),data,eps); % Find mols in neighborhood
        N = N(:); % N is an array of molecules in neighborhood     
        if numel(N) >= minp  % noise criteria
            C = C + 1; % If new increment cluster count
            clust(i) = C; % Apply cluster label
            
            hold on
            if C == 7
            plot(data(N,1),data(N,2),'.r')
            end
            while sum(clust(N) == 0) > 0 % while there are unvisited points
                ind = find(clust(N) == 0);
                ind2 = clust(N) == -1; % Incorporate all 'noise' variables into the cluster
                clust(ind2) = C; % Convert 'noise' to border points
                ID = N(ind(1)); % Look at first ungrouped cluster point
                
                ns = func_range_scan(data(ID,:),data,eps);
                if C == 7
                plot(data(ns,1),data(ns,2),'.g')
                plot(data(ID,1),data(ID,2),'.b')
                xlim([data(ID,1) - eps*5, data(ID,1) + eps*5])
                ylim([data(ID,2) - eps*5, data(ID,2) + eps*5])
                
                waitforbuttonpress
                end
                if numel(ns) >= minp
                    N = unique([N;ns]);
                    clust(ID) = C;
                    if C == 7
                    ind = find(clust(N) == 0);
                    plot(data(N,1),data(N,2),'.r')
                    plot(data(N(ind),1),data(N(ind),2),'.c')
                    waitforbuttonpress
                    end
                else
                    clust(ID) = -1;
                end
            end
        else % Minimum number of points not satified
            clust(i) = -1; % Set cluster to 'noise' cluster
        end
    end
    % movie2gif(M,'C:\Users\AJN Lab\Dropbox\Lab Meetings AJN\July 22 2019\DB_example.gif','DelayTime', 0.03,'LoopCount',Inf);
end


function ids = func_range_scan(me, them, eps)

dists = them(:,1)*0;
for j = 1:numel(them(1,:))
    dists = dists + (them(:,j) - me(j)).^2;
end
ids = find(dists <= eps^2);
end
end