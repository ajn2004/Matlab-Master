function [wnew, Cnew, clus_id] = web_cutter(w, C, data, dlims, R)
% This function will allow the user to select a distance to threshold
% connections with
wnew = w;
Cnew = C;

%% Data Representation
% Fill out this section if your data is representable as a graph

clust_loops(C,w, dlims(1),dlims(2),dlims(3), R)

clus_id = w(:,1).*0;
ask2 = 2000;
while true
    
    [row, col] = find(Cnew > 0);
    for i = 1:numel(row)
        dist = ((w(row(i),1)-w(col(i),1))^2 +(w(row(i),2)-w(col(i),2))^2)^0.5;
        if dist > ask2/1000
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
    
    %% Show Data
    subplot(2,2,2);plot(data(:,1),data(:,2),'.k'); % plot data
    hold on
    
    for i = 1:numel(Cnew(:,1))
        for j = i:numel(Cnew(1,:))
            if Cnew(i,j) == 1
                plot([w(i,1),w(j,1)],[w(i,2),w(j,2)],'b','LineWidth',1.5);
            end
        end
    end
    
    scatter(w(:,1),w(:,2),30,clus_id,'Filled'); % plot positions of nodes
    colormap('lines')
    
    hold off
    axis image
    drawnow
    
    disp(['There are ', num2str(clusts),' clusters found']);
    ask1 = input('Are you happy with this connection scheme? (y/n) :', 's');
    if strcmp(ask1,'y') || strcmp(ask1,'Y') % if the user is happy with the connection map
        break;
    end
    ask2 = input('What distance would you like to try? :');
    Cnew = C;
    [row, col] = find(Cnew > 0);
    
    for i = 1:numel(row)
        dist = ((w(row(i),1)-w(col(i),1))^2 +(w(row(i),2)-w(col(i),2))^2)^0.5;
        if dist > ask2/1000
            Cnew(row(i),col(i)) = 0;
            Cnew(col(i),row(i)) = 0;
        end
    end
    
    
    
    
    
    
    
    
end
end
