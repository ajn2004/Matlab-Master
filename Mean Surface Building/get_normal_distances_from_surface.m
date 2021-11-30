function distances = get_normal_distances_from_surface(centers, normals, data, cutoff_distance)


distances = zeros(numel(data(:,1)),2) + 20;
[Idx, D] = knnsearch(centers,data);

for i = 1:numel(data(:,1)) % loop over all red molecules
    % Find nearest normal center
    %     [Idx, D] = knnsearch(centers,[xfr(i),yfr(i),zfr(i)]);
    if D(i) < cutoff_distance % if nearest normal center is in range include measurement
        % build u and v vectors for distance calculation
        v = [data(i,1) - centers(Idx(i),1),data(i,2) - centers(Idx(i),2),data(i,3) - centers(Idx(i),3)];
%         u = [normals(Idx(i),1) + centers(Idx(i),1), normals(Idx(i),2) + centers(Idx(i),2),normals(Idx(i),3) + centers(Idx(i),3)];
        u = [normals(Idx(i),1), normals(Idx(i),2),normals(Idx(i),3)];
        % dot v w/ us to get overlap
        
        distances(i,1) = sum(v.*u);
%         if data(i,1) > 1.666 && data(i,1) < 1.668
%             a = 2;
%             b =5;
%         end
%         distances(i,1) = D(i);
        distances(i,2) = Idx(i);
    end
end
index = distances(:,1) < cutoff_distance;
% distances = distances(index,:);