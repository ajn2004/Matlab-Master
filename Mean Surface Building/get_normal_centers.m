function [centers, normals] = get_normal_centers(sdata,shrink)

k = boundary(sdata,shrink);
scatter3(sdata(:,1),sdata(:,2),sdata(:,3),[],'filled')
normals = [];
centers = [];
for i = 1:numel(k(:,1))
    v = [sdata(k(i,2),1) - sdata(k(i,1),1),sdata(k(i,2),2) - sdata(k(i,1),2),sdata(k(i,2),3) - sdata(k(i,1),3)];
    u = [sdata(k(i,3),1) - sdata(k(i,1),1),sdata(k(i,3),2) - sdata(k(i,1),2),sdata(k(i,3),3) - sdata(k(i,1),3)];
    normals(i,:) = cross(v,u); % find normal vector
    normals(i,:) = normals(i,:)/(sum(normals(i,:).*normals(i,:)))^0.5; % normalize to unit length
    centers(i,:) = mean(sdata(k(i,:),:)); % record central point
end