function [clust] = DB_scan(data, eps, minp)
% This function applies the DB scan algorithm

C = 1;
<<<<<<< HEAD
clust = (1:numel(data(:,1)))*0;

% scatter(data(:,1),data(:,2),2,[0 0 0],'Filled')
count = 1;
for i = 1:numel(data(:,1))
    
    if clust(i) == 0 % If we haven't visited the point yet, perform analysis
        N = func_range_scan(data(i,:),data,eps); % Find mols in neighborhood
        N = N(:); % N is an array of molecules in neighborhood
=======
visit = (1:numel(data(:,1)))*0;
clust = visit;
scatter(data(:,1),data(:,2),2,[0 0 0],'Filled')
count = 1;
for i = 1:numel(data(:,1))
    if visit(i) == 0 % If we haven't visited the point yet, perform analysis
        visit(i) = 1; % mark as visited
        N = func_range_scan(data(i,:),data,eps); % Find mols in neighborhood
        N = N(:);
>>>>>>> master

        if numel(N) >= minp  % noise criteria
            clust(i) = C;
            C = C + 1;
<<<<<<< HEAD

        
            test = clust(N) == 0; % see if anything in neighborhood is unvisited
            while sum(test) > 0 % while there are unvisited points
                ind = find(clust(N) == 0);
                ind2 = clust(N) == -1; % Incorporate all 'noise' variables into the cluster
                clust(ind2) = C;
                ID = N(ind(1));
                
                ns = func_range_scan(data(ID,:),data,eps);
                if  clust(ID) == 0 % visit criteria
                    if numel(ns) >= minp % noise criteria
                        N = unique([N;ns(:)]);
                        clust(ID) = C; % mark point as apart of cluster
                    else
                        clust(ID) = -1; % mark point as apart of noise
                    end
                end
                
                test = clust(N) == 0;
            end
        else % Minimum number of points not satified
            clust(i) = -1; % Set cluster to 'noise' cluster
        end
    end
end
% movie2gif(M,'C:\Users\AJN Lab\Dropbox\Lab Meetings AJN\July 22 2019\DB_example.gif','DelayTime', 0.03,'LoopCount',Inf);
end

function ids = func_range_scan(me, them, eps) 
% Algorithm to find molecules in neighborhood
=======
                    hold on
        scatter(data(i,1),data(i,2),5,clust(i),'Filled');
        scatter(data(N,1),data(N,2),3,ones(numel(N),1)*clust(i),'Filled');
        
            test = visit(N) == 0; % see if anything in neighborhood is unvisited
            while sum(test) > 0 % while there are unvisited points
                ind = find(visit(N) == 0);
                ID = N(ind(1));
                visit(ID) = 1; % mark point as visited
                ns = func_range_scan(data(ID,:),data,eps);
                scatter(data(ID,1),data(ID,2),3,[1 0 0])
                scatter(data(ns,1),data(ns,2),3,ones(numel(ns),1)*clust(i),'Filled');
                M(count) = getframe(gcf);
                count = count +1;
                drawnow
                if numel(ns) >= minp  % noise criteria
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
movie2gif(M,'C:\Users\AJN Lab\Dropbox\Lab Meetings AJN\July 22 2019\DB_example.gif','DelayTime', 0.03,'LoopCount',Inf);
end

function ids = func_range_scan(me, them, eps)
>>>>>>> master
dists = them(:,1)*0;
for i = 1:numel(them(1,:))
    dists = dists + (them(:,i) - me(i)).^2;
end
ids = find(dists <= eps^2);
end