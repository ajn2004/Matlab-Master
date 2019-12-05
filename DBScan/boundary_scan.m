function [cluster, clust] = boundary_scan(data, eps, minp, scal)
% This function applies the DB scan algorithm and follows it up by a matlab
% boundary 

C = 1;
visit = (1:numel(data(:,1)))*0;
clust = visit;
for i = 1:numel(data(:,1))
    if visit(i) == 0 % If we haven't visited the point yet, perform analysis
        visit(i) = 1; % mark as visited
        N = func_range_scan(data(i,:),data,eps); % Find mols in neighborhood
        N = N(:);
        if numel(N) >= minpts  % noise criteria
            clust(i) = C;
            C = C + 1;
            
            test = visit(N) == 0; % see if anything in neighborhood is unvisited
            while sum(test) > 0 % while there are unvisited points
                ind = find(visit(N) == 0);
                ID = N(ind(1));
                visit(ID) = 1; % mark point as visited
                ns = func_range_scan(data(ID,:),data,eps);
                if numel(ns) >= minpts  % noise criteria
                    N = [N;ns(:)];
                    
                end
                if clust(ID) == 0
                    clust(ID) = clust(i);
                end
                test = visit(N) == 0;
            end
        end
    end
end

for i = 1:max(clust)
    ind = find(clust == i);
    bound = boundary(xf(ind),yf(ind),zf(ind));
    cluster(i).IDs = ind(bound);
end
end

function ids = func_range_scan(me, them, eps)
dists = them(:,1)*0;
for i = 1:numel(them(1,:))
    dists = dists + (them(:,i) - me(i)).^2;
end
ids = find(dists <= eps^2);
end