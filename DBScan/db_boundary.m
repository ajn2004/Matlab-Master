%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary Definer
% This algorithm allows you to select a file to perform a DB scan followed
% by a boundary identifcation for futher study of structures in superres
% data sets
%
% AJN 6/26/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
clearvars;

eps = 0.2; % minimum distance in microns
min_pts = 4; %
scl = 0.75;
r = 0;

[fname, fpath] = uigetfile('*mat');
load([fpath,fname],'xf_fixed','yf_fixed','ncoords','q','framenumber');
% Setup Position
xf = xf_fixed*q;
yf = yf_fixed*q;
zf = func_shift_correct(ncoords(:,3)*q,framenumber,r);
data = [xf,yf,zf];

clust = DB_scan(data,eps,min_pts);
for i = 1:numel(clust)
    ind1 = find(clust == i);
    cluster(i).bound = ind1(boundary(xf(ind1),yf(ind1),zf(ind1),scl));
    
    if ~isempty(cluster(i).bound)
        for j = 1:numel(cluster(i).bound(:,1))
            ind = [cluster(i).bound(j,:),cluster(i).bound(j,1)];
            plot3(xf(ind),yf(ind),zf(ind),'k')
            hold on
        end
        plot3(xf(ind1),yf(ind1),zf(ind1),'.');
    end
end
hold off
axis equal