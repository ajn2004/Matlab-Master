%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 Color transform
% 
% My attempt at making a 2-Image Region correlator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all
clc;

pix_thresh = 40;

files = dir('*dast*');
xf = [];
yf = [];
N = [];
fnum = [];
for i = 1:numel(files)
    load(files(i).name)
    xf = [xf; fits(:,1)];
    yf = [yf; fits(:,2)];
    N = [N; fits(:,3)];
    fnum = [fnum; framenumber + max([max(fnum),0])];
end
id = N < 200;
xf(id) = [];
yf(id) = [];
fnum(id) = [];
plot3(xf,yf,fnum,'.');

split = 180;
% hold on
% plot(split,6,'.r')
% legend('Data','Split Point')

% Split localizations based on position
spind = xf < split;
xfb = xf(spind);
xfr = xf(logical(1-spind));
yfb = yf(spind);
yfr = yf(logical(1-spind));
fnumb = fnum(spind);
fnumr = fnum(logical(1-spind));
plot(xfb,yfb,'b.')
hold on
plot(xfr,yfr,'r.');
hold off
title('Channel Identification')
figure
xnb = xfb - mean(xfb);
ynb = yfb - mean(yfb);
xnr = xfr - mean(xfr);
ynr = yfr - mean(yfr);
mfa = xnr.*0;
plot(xnb,ynb,'.b')
hold on
plot(xnr,ynr,'.r')
hold off
axis equal
title('Overlay')
figure
pairs = [yfb]*0; % pairs will link blue to red, so the index of pairs will correspond to blue and the value at that index corresponds to the indexed red molecule

% pair molecules based on frame and nearest distance
% Start by looking at 1 molecule in a frame, then filter by frame number
% Then filter by distance
for i = 1:numel(xnb)
    fids = find(fnumr == fnumb(i)); % Find all ID's w/ same frame number   
    dists = ((xnb(i) - xnr(fids)).^2 + (ynb(i) - ynr(fids)).^2).^0.5;
    id = find(dists == min(dists));
    if min(dists) <= pix_thresh
        id1 = find(pairs == fids(id));
        if numel(id1) >= 1 % if red molecule is already paired, only accept closest distance pair
            odist = ((xnb(id1) - xnr(fids(id))).^2 + (ynb(id1) - ynr(fids(id))).^2).^0.5;
            if odist > min(dists)
                pairs(id1) = -1;
                pairs(i,1) = fids(id);
                mfa(fids(id)) = i;
            else
                pairs(i,1) = -1;
            end
        else % there were no other pairs
            pairs(i,1) = fids(id);
            mfa(fids(id)) = i;
        end
    else
        pairs(i,1) = -1;
    end
end
deleteid = find(pairs == -1);

xnb(deleteid) = [];
xfb(deleteid) = [];
ynb(deleteid) = [];
yfb(deleteid) = [];
pairs(deleteid) = [];
% for i = 1:numel(pairs)
%     ind = find(pairs == pairs(i));
%     if numel(ind) > 1
%         dists = ((xnb(ind) - xnr(pairs(i))).^2 + (ynb(i) - ynr(fids)).^2).^0.5;
%         id = dists ~= min(dists);
%         pairs(ind(id)) = [];
%         xnb(ind(id)) = [];
%         xfb(ind(id)) = [];
%         xnb(ind(id)) = [];
%     end
% end



% for i = 1:max(fnum)
%     indb = find(fnumb == i);
%     indr = find(fnumr == i);
%     for j = 1:numel(indb)
%         dists = ((xnb(indb(j)) - xnr(indr)).^2 + (ynb(indb(j)) - ynr(indr)).^2).^0.5;
%         id = find(dists == min(dists));
%         pairs(indb(j),1) = indr(id);
%     end
% end


xpr = xnr(pairs);
ypr = ynr(pairs);
% mag = mean((xpr.^2 + ypr.^2).^0.5./(xnb.^2 + ynb.^2).^0.5);
ins = [xfb, yfb];
outs = [xpr,ypr];

inds = randperm(numel(outs(:,1)));

xsub = ins;
% x = ins(inds(round(0.7*numel(inds))):end,:);
ysub =  outs;


xnet = feedforwardnet([50]);
% xnet.Inputs{1}.size = 2;
ynet = xnet;
xnet = train(xnet,xsub.',ysub.');
% xnet = train(xnet,xsub.',ysub(:,1).');
% ynet = train(ynet,xsub.',ysub(:,2).');
xs = xnet(xsub.');
x0 = xs(1,:);
y0 = xs(2,:);
xo = x0(:);
yo = y0(:);
% xo = xnet(xsub.');
% yo = ynet(xsub.');
xo = xo(:);
yo = yo(:);
% save('Channel_net.mat','xnet','ynet');
plot(xo,yo,'.b')
hold on
plot(xnr,ynr,'.r')
legend('Transformed','2nd Channel')
