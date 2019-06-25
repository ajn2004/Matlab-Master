%% DB SCAN
% A standalone program to build the DB SCan algorithm for 3D localization
% Data
% AJN 6/18/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;

% Grab a test set
fname = 'Cell2_dz10_r1_1_dast_tol_dc_150nm_traj.mat';
fpath = 'C:\Users\AJN Lab\Dropbox\Data\6-3-19 cono-halo Neurons\Analysis\toleranced\DC\traj\';
file = [fpath, fname];
load(file,'ncoords','xf_fixed','yf_fixed','q','framenumber');

eps = 0.2;
minpts = 4;

% Setup Position
xf = xf_fixed*q;
yf = yf_fixed*q;
zf = func_shift_correct(ncoords(:,3)*q,framenumber,1).';
C = 1;
visit = (1:numel(xf))*0;
clust = visit;
for i = 1:numel(xf)
    if visit(i) == 0 % If we haven't visited the point yet, perform analysis
        visit(i) = 1; % mark as visited
        N = func_range_scan([xf(i),yf(i),zf(i)],[xf,yf,zf],eps); % Find mols in neighborhood
        N = N(:);
        if numel(N) >= minpts  % noise criteria
            clust(i) = C;
            C = C + 1;
            
            test = visit(N) == 0; % see if anything in neighborhood is unvisited
            while sum(test) > 0 % while there are unvisited points
                ind = find(visit(N) == 0);
                ID = N(ind(1));
                visit(ID) = 1; % mark point as visited
                ns = func_range_scan([xf(ID),yf(ID),zf(ID)],[xf,yf,zf],eps);
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

plot3(xf,yf,zf,'.k','MarkerSize',4)
hold on

for i = 1:max(clust)
    ind = find(clust == i);
    plot3(xf(ind),yf(ind),zf(ind),'.','MarkerSize',10);
    bound = boundary(xf(ind),yf(ind),zf(ind));
    cluster(i).bound = ind(bound);
end

hold off
axis equal
figure
for i = 1:numel(cluster)
    if ~isempty(cluster(i).bound)
    for j = 1:numel(cluster(i).bound(:,1))
        ind = [cluster(i).bound(j,:),cluster(i).bound(j,1)];
        plot3(xf(ind),yf(ind),zf(ind),'k')
        hold on
    end
    end
end
hold off 
axis equal

function ids = func_range_scan(me, them, eps)
dists = them(:,1)*0;
for i = 1:numel(them(1,:))
    dists = dists + (them(:,i) - me(i)).^2;
end
ids = find(dists <= eps^2);
end