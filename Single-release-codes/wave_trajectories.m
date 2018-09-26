% Wave Trajectories
% Display trajectories of wavelet analysis localizations
% close all
% xf = fits(:,1) + cents(:,1);
% yf = fits(:,2) + cents(:,2);
% N = fits(:,3);

dmax = 300;
% ind  = fits(:,3) > 70;
% ind  = N > 200;
xf_all = xf(ind);
yf_all = yf(ind);
framenum_all = fnum(ind);
trajec = struct('t',{[]}); % initialize trajectory variable
foll = zeros(numel(framenum_all),1);
% xc = fits(ind,1);
% yc = fits(ind,2);
% sdi1 = iloc(:,:,ind);
% iloc = sdi1;
%% This has been copied from 'traj_analysis.m'
% loop over all frames to build the connections
for i = 1:max(framenum_all)
    clear dist
    % i is the framenumber, so we want to look at
    if mod(i,2) == 1 % special conditions for frame 1 because there is no previous
        cind = find(framenum_all == i); % current index for loop i
        fodex = find(framenum_all == i+1); % index of all molecules on following frame
        if ~isempty(cind) && ~isempty(fodex)
            for j = 1:numel(cind)
                for k = 1:numel(fodex)
                    dist(k) = q*1000*((xf_all(cind(j)) - xf_all(fodex(k))).^2 + (yf_all(cind(j)) - yf_all(fodex(k))).^2).^0.5;
                end
                disto = q*1000*((xf_all(cind(j)) - xf_all(cind)).^2 + (yf_all(cind(j)) - yf_all(cind)).^2).^0.5; 
                flag = find(disto > 0 & disto <= 2*dmax);
                
                if min(dist) <= dmax && sum(dist <= 2*dmax) < 2 && isempty(flag)
                    foll(cind(j)) = fodex(find(dist == min(dist)));
                end
            end
        end
    end
end

%% loop over conections to build trajectory variable
count = 1;
for i = 1:numel(foll) % foll is an array of follows
    traj = []; % initiallize traj number
    if foll(i) ~= 0 % if there is a following index
        traj = i; % Start traj w/ current index
        [traj, foll] = traj_connect(foll,traj, foll(i)); % send all follows, the traj you're building, and the current follower
    trajec(count) = struct('t', traj);
  
    count = count+1;
    end
    
end

set_scale(std(dip1,1,3),q,4);
colormap('gray');
hold on
ind1 = mod(framenum_all,2) == 1;


plot(xf_all(ind1),yf_all(ind1),'.g')
plot(xf_all(logical(1-ind1)),yf_all(logical(1-ind1)),'.r')
for i = 1:numel(trajec)    
    inds = trajec(i).t;
    plot(xf_all(inds),yf_all(inds),'b');
end
legend('Releases','Post-stims','Trajectory')
hold off
axis image

% % show all localizations that got through the threshold
% f = figure;
% tg = uitabgroup(f);
% t5 = uitab(tg,'Title','Frames');
% t55 = uitabgroup(t5);
% o = numel(xc);
% [m,n] = size(sdi1(:,:,1));
% for i = 1:o
%     ax = axes(uitab(t55,'Title',['F ', num2str(i)]));
%     imagesc(ax,sdi1(:,:,i));
%     hold on
% %     wind = -3:3;
%     colormap(ax,'gray')
% %     [xf,yf,sx,sy,Neat,O] = mle_Gauss(sdi1((m+1)/2+wind,(n+1)/2+wind,i),110*pi/180);
% %     xfm(i) = xf + cents(i,1);
% %     yfm(i) = yf + cents(i,2);
% %         plot(ax,xf_all(i) - cents(i,1),yf_all(i)- cents(i,2),'bx');
%     plot(ax,xc(i)+(n+1)/2,yc(i)+ (m+1)/2,'rx');
% %     plot(xc(i)+,yc(i),'r.')
% %     title(['Sx = ', num2str(sx),' Sy = ', num2str(sy), 'N = ', num2str(Neat),' O = ', num2str(O)]);
%     hold off
%     axis image
% end