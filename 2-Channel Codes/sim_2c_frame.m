% Sim 2 channel Frame
close all
clearvars
clc

tx = 12; % horizontal translation in um
ty = -0.5; % vertical translation in um
rz = 2; % image rotation between frames in degrees
 m = 1.05; % magnification between images
 
 n = 100000; % number of molecule pairs to consider
 fms = 200000;
 
 lp = 5;  % nm precision of localization
 
 theta = deg2rad(rz);
 % build localizations
for i = 1:n
    xt = 10*rand(1)+1;
    yt = 10*rand(1)+1;
    xf1(i,1) = xt + lp/1000*randn(1);
    yf1(i,1) = yt + lp/1000*randn(1);
    tpair(i,1) = i;
    xf2(i,1) = m*(xt*cos(theta) - yt*sin(theta) + lp/1000*randn(1)) + tx;
    yf2(i,1) = m*(yt*cos(theta) + xt*sin(theta) + lp/1000*randn(1)) + ty;
    fnum(i,1) = floor(i*fms/n)+1;
end

plot(xf1,yf1,'.b')
hold on
plot(xf2,yf2,'.r')
title('Simulated Data')
legend('Blue Channel','Red Channel')
hold off
axis equal

% Shuffle localizations
xl = [xf1;xf2];
yl = [yf1;yf2];
tps = [tpair;tpair];
fn = [fnum;fnum];
dlist = randperm(numel(xl));
xf = xl(dlist);
yf = yl(dlist);
fnum = fn(dlist);
truths = tps(dlist);

% figure
% plot(xf,yf,'.')
% axis equal
% title('Unidentified locs')

split = 11.5;
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
trb = truths(spind);
trr = truths(logical(1-spind));
figure
% plot(xfb,yfb,'.b');
% axis equal
% hold on
% plot(xfr,yfr,'.r');
% hold off
% title('Data after Split')
% center localizations based on mean
xnb = xfb - mean(xfb);
ynb = yfb - mean(yfb);
xnr = xfr - mean(xfr);
ynr = yfr - mean(yfr);
% plot(xnb,ynb,'.b')
% hold on
% plot(xnr,ynr,'.r')
% hold off
% axis equal

pairs = [yfb]*0; % pairs will link blue to red, so the index of pairs will correspond to blue and the value at that index corresponds to the indexed red molecule

% pair molecules based on frame and nearest distance
for i = 1:fms+1
    indb = find(fnumb == i);
    indr = find(fnumr == i);
    for j = 1:numel(indb)
        dists = ((xnb(indb(j)) - xnr(indr)).^2 + (ynb(indb(j)) - ynr(indr)).^2).^0.5;
        id = find(dists == min(dists));
        pairs(indb(j),1) = indr(id);
    end
end
% sz= 4;
% scatter(xnb,ynb,sz,pairs);
% hold on
% scatter(xnr(pairs),ynr(pairs),sz,pairs);
% hold off
% colormap('hsv')
% axis equal
% title('Overlay of two channels');
tgr = trr(pairs);

% figure
% for i = 1:360
%     ang = deg2rad(i);
%     scatter(xnb,ynb,sz,pairs);
%     hold on
%     scatter(xnr(pairs)*cos(ang) - ynr(pairs)*sin(ang), xnr(pairs)*sin(ang) + cos(ang)*ynr(pairs),sz,pairs)
%     colormap('jet')
%     xlim([-7 7])
%     ylim([-7 7])
%     M(i) = getframe(gcf);
%     hold off
% end
sum(tgr ~= trb)

xpr = xnr(pairs);
ypr = ynr(pairs);
% mag = mean((xpr.^2 + ypr.^2).^0.5./(xnb.^2 + ynb.^2).^0.5);
ins = [xnb, ynb];
outs = [xpr,ypr];

inds = randperm(numel(outs(:,1)));

xsub = ins(inds(1:round(0.7*numel(inds))),:);
x = ins(inds(round(0.7*numel(inds))):end,:);
ysub =  outs(inds(1:round(0.7*numel(inds))),:);
tsts = outs(inds(round(0.7*numel(inds))):end,:);

xnet = feedforwardnet([20]);
xnet.Inputs{1}.size = 2;
ynet = xnet;

xnet = train(xnet,xsub.',ysub(:,1).');
ynet = train(ynet,xsub.',ysub(:,2).');

xo = xnet(x.');
yo = ynet(x.');
xo = xo(:);
yo = yo(:);

plot(x(:,1),x(:,2),'.')
hold on
plot(xo,yo,'.');
hold off
figure
histogram(tsts(:,1)-xo)
hold on
histogram(tsts(:,2) - yo)
hold off
xlabel(['X uncertainty is', num2str(std(tsts(:,1)-xo)*1000), 'nm'])
ylabel(['Y uncertainty is', num2str(std(tsts(:,2)-yo)*1000), 'nm'])
title(['Simulated Precision is ', num2str(lp),'nm'])
% [X, Y] = meshgrid(-5:0.1:5);
% lx = X(:);
% ly = Y(:);
% xl = xnet([lx,ly].');
% yl = ynet([lx,ly].');
% xl = xl(:);
% yl = yl(:);
% plot(lx,ly,'.')
% hold on
% plot(xl,yl,'.')