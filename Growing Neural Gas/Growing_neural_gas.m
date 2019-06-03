% Growing Neural Gas
% Based off the algorithm described by Fritzke "A Growing Neural Gas
% Network Learns Topologies"
% AJN 4-9-19
clearvars;
close all;
clc;
%
epsb = 0.2;
epsn = 0.006;
amax = 50;
lamb = 100;
alpha = 0.5;
decr = 0.995;
vangle = [00, 90];
% [fname, fpath] = uigetfile('*mat');
file = ['C:\Users\AJN Lab\Dropbox\Data\5-22-19 Storm Neurons\Analysis\toleranced\DC\traj\Cell3_dz0_r0_10mw_dast_tol_dc_150nm_traj.mat'];
fname = 'hek6_r2_dz20_dast_tol_dc_100nm_traj.mat';
load(file,'ncoords','framenumber','q');

shows = [1:20:1000,1000:200:10000,10000:2000:50000];
% shows = [];
% Localization based preallocation
r = str2num(fname(strfind(fname,'_r')+2));
zf = func_shift_correct(ncoords(:,3)*q,framenumber,r).';
xf = q*ncoords(:,1);
yf = q*ncoords(:,2);
N = numel(xf);

% Defining data vector
data = [xf,yf,zf];

% Initializing learning vector Step 0
units = 2;

T = [0 1; 1 0];
A = [-1 0; 0 -1];
id = randi(N);
w = [xf(id),yf(id),zf(id)];
id = randi(N);
w = [w;[xf(id),yf(id),zf(id)]];
del = [0,0];
count = 0;
% plot3(data(:,1),data(:,2),data(:,3),'.k')
figure('Units','Normalized','OuterPosition',[0 0 1 1])
while true
    % Step 1 generate random signal
    id = randi(N);
    count = count +1;
    % Step 2 Find nearest and 2nd nearest units
    d = 0;
    dist = [];
    for i = 1:units
        d = 0;
        for j = 1:numel(data(1,:))
            d = d + (w(i,j) - data(id,j))^2;
        end
        d = d.^0.5;
        dist(i) = d;
    end
    
    % Step 3 Increment age of all connections to nearest node by 1
    [rank_dist, backrank, ranking] = unique(dist);
    k = ranking - 1;
    io = find(k == 0);
    i1 = find(k == 1);
    A(io,:) = A(io,:) + 1.*(A(io,:) > -1);
    A(:,io) = A(:,io) + 1.*(A(:,io) > -1);
    
    % Step 4 Increment Error between nearest point and signal point
    d = 0;
    for j = 1:numel(data(1,:))
        d = d + (w(io,j) - data(id,j))^2;
    end
    del(io) = del(io) + d;
    
    % Step 5 Move closest and neighbors towards signal by set fractions
    for i = 1:numel(data(1,:))
        w(io,i) = w(io,i) + epsb*(data(id,i) - w(io,i));
    end
    
    ind = find(T(io,:) == 1);
    for j = ind
        for i = 1:numel(data(1,:))
            w(j,i) = w(j,i) + epsn*(data(id,i) - w(j,i));
        end
    end
    
    % Step 6 Connect / age nearest and 2nd nearest points by an edge
    if T(io,i1) == 1
        A(io,i1) = 0;
        A(i1,io) = 0;
    else
        T(io,i1) = 1;
        T(i1,io) = 1;
        A(io,i1) = 0;
        A(i1,io) = 0;
    end
    
    %Step 7 Remove connections older than amax
    [row,col] = find(A >= amax);
    for i = 1:numel(row)
        T(row(i),col(i)) = 0;
        T(col(i),row(i)) = 0;
        A(row(i),col(i)) = -1;
        T(col(i),row(i)) = -1;
    end
    
    % Find nodes w/out connections and remove them from the set
    st = sum(T,1);
    inds = find(st == 0);
    if ~isempty(inds)
        T(inds,:) = [];
        T(:,inds) = [];
        A(inds,:) = [];
        A(:,inds) = [];
        w(inds,:) = [];
        del(inds) = [];
        units = units - numel(inds);
    end
    
    % Step 8
    if mod(count,lamb) == 0 % if we are an integer multiple of lambda
        count
        [rank_dist, backrank, ranking] = unique(del);
        k = find(ranking == max(ranking),1);
        % K constitutes largest error member of the set
        % Find K's neighbor who has largest error
        inds = find(T(k,:) == 1);
        indy = find(del(inds) == max(del(inds)),1);
        i1 = inds(indy);
        
        % Insert new unit between these two neighbors
        for j = 1:numel(data(1,:))
            w(units + 1,j) = 0.5*(w(k,j) + w(i1,j));
        end
        
        % Update connections appropriately
        units = units +1;
        T(k,i1) = 0;
        T(i1,k) = 0;
        A(k,i1) = -1;
        A(i1,k) = -1;
        T(k,units) = 1;
        T(units,k) = 1;
        T(i1,units) = 1;
        T(units,i1) = 1;
        A(k,units) = 0;
        A(units,k) = 0;
        A(:,units) = -1;
        A(units,:) = -1;
        A(k,units) = 0;
        A(units,k) = 0;
        A(units,i1) = 0;
        A(i1,units) = 0;
        del(k) = alpha*del(k);
        del(i1) = alpha*del(i1);
        del(units) = del(k);
        %     % Visualization
        %     plot3(data(:,1),data(:,2),data(:,3),'.')
        %     hold on
        %     plot3(w(:,1),w(:,2),w(:,3),'.r')
        %     axis equal
        %     for i = 1:units
        %         ind = find(T(i,:) == 1);
        %             for j = 1:numel(ind)
        %                 plot3([w(i,1),w(ind(j),1)],[w(i,2),w(ind(j),2)],[w(i,3),w(ind(j),3)],'b')
        %             end
        %     end
        %     drawnow
        %     hold off
    end
    % Step 9 Decrease all error by constant d
    del = decr*del;
    
    if ismember(count,shows)
%         shows = shows + 1;
        subplot(1,3,1)
        plot3(data(:,1),data(:,2),data(:,3),'k.','MarkerSize',1)
        view(vangle)
        axis equal
        title('Raw Data')
        subplot(1,3,2);
        plot3(data(:,1),data(:,2),data(:,3),'k.','MarkerSize',3)
        
        %         axis equal
        hold on
        plot3(w(:,1),w(:,2),w(:,3),'.r','MarkerSize',10)
        %         view([-183,38])
        view(vangle)
        axis equal
        for i = 1:units
            ind = find(T(i,:) == 1);
            for j = 1:numel(ind)
%                 plot3([w(i,1),w(ind(j),1)],[w(i,2),w(ind(j),2)],[w(i,3),w(ind(j),3)],'y','LineWidth',3)
                plot3([w(i,1),w(ind(j),1)],[w(i,2),w(ind(j),2)],[w(i,3),w(ind(j),3)],'b','LineWidth',2)
            end
        end

        hold off
        title('Data w/ Neural Gas')
        
        subplot(1,3,3)
        plot3(w(:,1),w(:,2),w(:,3),'.r','MarkerSize',10)
        view(vangle)
        hold on
        axis equal
        title('Neural Gas with Connections')
        title(['Iteration ', num2str(count),])
        for i = 1:units
            ind = find(T(i,:) == 1);
            for j = 1:numel(ind)
                plot3([w(i,1),w(ind(j),1)],[w(i,2),w(ind(j),2)],[w(i,3),w(ind(j),3)],'y','LineWidth',3)
                plot3([w(i,1),w(ind(j),1)],[w(i,2),w(ind(j),2)],[w(i,3),w(ind(j),3)],'b','LineWidth',1)
            end
        end
        hold off
                drawnow
        if exist('M') ~= 1
            M(1) = getframe(gcf);
        else
            M(numel(M)+1) = getframe(gcf);
        end
    end
    if count > 50000
        break
    end
end
% movie2gif(M,'Results_of_GNG.gif','DelayTime',0.1);
% for i = 1:360
%     view([-183+i,38])
%     drawnow
%     M(i) = getframe(gcf);
% end
% movie2gif(M,'Final_gng.gif','DelayTime',0.1);