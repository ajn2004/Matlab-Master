function surface_data = get_surface_from_data(data, shrink, max_dist)
k = boundary(data,shrink);
u_k = unique(k);
smooth_data = [];
% Grab neighbors of outter-most localizations to incorporate in surface
% calculation
surface_data = [];
for i = 1:numel(u_k)
    % Loop over all surface molecules
    %     index = ((xfr - xfr(u_k(i))).^2 + (yfr - yfr(u_k(i))).^2 + (zfr - zfr(u_k(i))).^2 ).^0.5 <= max_dist; % find molecules within range
    index = ((data(:,1) - data(u_k(i),1)).^2 + (data(:,2) - data(u_k(i),2)).^2 + (data(:,3) - data(u_k(i),3)).^2 ).^0.5 <= max_dist; % find molecules within range
    surface_data = [surface_data;[data(index,1),data(index,2),data(index,3)]]; % add in range molecules to surface data
end
surface_data = unique(surface_data,'rows'); % delete duplicates