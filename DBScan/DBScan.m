%% DB SCAN
% A standalone program to build the DB SCan algorithm for 3D localization
% Data
% AJN 6/18/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;

% Grab a test set
[fname, fpath] = uigetfile('*mat');
file = [fpath, fname];
load(file,'ncoords','xf_fixed','yf_fixed','q','framenumber');

eps = 0.225;
minpts = 7;
scle = 0.8;

% Setup Position
xf = xf_fixed*q;
yf = yf_fixed*q;
zf = func_shift_correct(ncoords(:,3)*q,framenumber,1);
C = 1;
visit = (1:numel(xf))*0;
clust = visit;
% plot(xf,yf,'.k');
% hold on
for i = 1:numel(xf)
    if visit(i) == 0 % If we haven't visited the point yet, perform analysis
        visit(i) = 1; % mark as visited
        N = func_range_scan([xf(i),yf(i),zf(i)],[xf,yf,zf],eps); % Find mols in neighborhood
        N = N(:);
%         plot(xf(i),yf(i),'.r','MarkerSize',5)
%         plot(xf(N),yf(N),'.g','MarkerSize',3);
%         drawnow
        if numel(N) >= minpts  % noise criteria
            clust(i) = C;
            C = C + 1;
            
            test = visit(N) == 0; % see if anything in neighborhood is unvisited
            while sum(test) > 0 % while there are unvisited points
                ind = find(visit(N) == 0);
                ID = N(ind(1));
                visit(ID) = 1; % mark point as visited
                ns = func_range_scan([xf(ID),yf(ID),zf(ID)],[xf,yf,zf],eps);
%                 plot(xf(ns),yf(ns),'.g','MarkerSize',3);
%                 drawnow
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

scatter3(xf,yf,zf,ones(numel(zf),1)*4,clust,'Filled')
hold on

for i = 1:max(clust)
    ind = find(clust == i);
    bound = boundary(xf(ind),yf(ind),zf(ind),scle);
    cluster(i).bound = ind(bound);
end

hold off
axis equal
figure
for i = 1:numel(cluster)
    if ~isempty(cluster(i).bound) %&& i == 2 || i == 8
    for j = 1:numel(cluster(i).bound(:,1))
        ind = [cluster(i).bound(j,:),cluster(i).bound(j,1)];
        plot3(xf(ind),yf(ind),zf(ind),'k')
        hold on
    end
    end
end
 
scatter3(xf,yf,zf,ones(numel(zf),1)*4,clust,'Filled')
axis equal
hold off
xlabel('Microns')
ylabel('Microns')
zlabel('Microns')
function ids = func_range_scan(me, them, eps)
dists = them(:,1)*0;
for i = 1:numel(them(1,:))
    dists = dists + (them(:,i) - me(i)).^2;
end
ids = find(dists <= eps^2);
end