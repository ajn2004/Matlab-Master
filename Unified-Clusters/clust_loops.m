function clust_loops(C,w, ds,dr,de,R)

count = 1;
for d = ds:dr:de
    d
    Cnew = C;
    [row, col] = find(Cnew > 0);
    
    for i = 1:numel(row)
        dist = ((w(row(i),1)-w(col(i),1))^2 +(w(row(i),2)-w(col(i),2))^2)^0.5;
        if dist > d/1000
            Cnew(row(i),col(i)) = 0;
            Cnew(col(i),row(i)) = 0;
        end
    end
      
    % Cluster Analysis
    clus_id = w(:,1).*0;
    clusts = 0;
    for i = 1:numel(Cnew(:,1)) % loop through rows
        if clus_id(i) == 0; % if a row is 0 add 1 to the clusts and assign to that cluster id
            clusts = clusts +1;
            clus_id(i) = clusts;
            for j = i:numel(Cnew(1,:)) % loop over column indecies to find where links are
                if Cnew(i,j) == 1 && clus_id(j) == 0
                    clus_id(j) = clus_id(i);
                    [Cnew, clus_id] = changerows(Cnew, j, clus_id); % run algorithm to go change rows of found indices
                end
            end
        end

    end
    clusters(count) = clusts;
    count = count +1;
end

subplot(2,2,3); plot(R, clusters); 
title('Number of clusters versus cutoff linkage');
xlabel('Distance');
ylabel('Clusters');
subplot(2,2,4); plot(diff(R)+R(1:end-1),diff(clusters));
title('Differential of number of clusters versus cutoff linkage');
xlabel('Distance');
ylabel('dClusters');