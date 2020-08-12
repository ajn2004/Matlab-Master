function clus_coords = get_z_tilt_parameters(data, points)
% Expecting x-y-z data in the form of [x,y,z] need to return a model of the
% data that can be used as a spline interpolation later on

% kmeans
[clus_coords] = funkmeanscluster(points,data);
% plot3(clus_coords(:,1),clus_coords(:,2),clus_coords(:,3),'.r','MarkerSize', 20)
% hold on
% plot3(xt,yt,zt,'.')
% hold off
