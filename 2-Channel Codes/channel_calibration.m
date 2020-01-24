clearvars
files = dir('*dast*'); % List calibration data files
% allocate arrays for building
xf = [];
yf = [];
N = [];
fnum = [];
% load calibration data and store in conglomorate array
for i = 1:numel(files)
    load(files(i).name)
    xf = [xf; fits(:,1)];
    yf = [yf; fits(:,2)];
    N = [N; fits(:,3)];
    fnum = [fnum; framenumber + max([max(fnum),0])];
end

% display calibration parameters of splitting channels via horizontal
% barrier
plot(xf,yf,'.')
split = 180;
hold on
plot([split,split],[min(yf),max(yf)],'r')
added = xf*0;

% Create one point for each set of points
X = [];
Y = [];
stdx = X;
stdy = Y;
% average points based on several frames of the same latex bead
for i = 1:numel(xf)
    if added(i) == 0
        dist = ((xf-xf(i)).^2 + (yf-yf(i)).^2 ).^0.5;
        ind = dist < 1;
        X = [X;mean(xf(ind))];
        Y = [Y;mean(yf(ind))];
        stdx = [stdx, std(xf(ind))];
        stdy = [stdy, std(yf(ind))];
        added(ind) = 1;
    end
end

% display averages and split data arrays
plot(X,Y,'.g')
id = X < split;
XR = X(id);
YR = Y(id);
XO = X(logical(1-id));
YO = Y(logical(1-id));
id = YR < min(YO);
XR(id) = [];
YR(id) = [];
close all
plot(XR,YR,'.r');
hold on
plot([split,split],[min(yf),max(yf)],'k');
plot(XO,YO,'.b')
% Idea of playing with boundary as a marker of the transformation space
bound_red = boundary(XR,YR);
bound_orange = boundary(XO,YO);
ind = [1 4 7 9 10 11 2 3 12 5];
XT = XO(bound_orange(ind));
YT = YO(bound_orange(ind));
XRT =XR(bound_red(ind));
YRT =YR(bound_red(ind));

plot(XO(bound_red),YO(bound_red),'r')
MX = xy_feature(XO(bound_orange(ind)),YO(bound_orange(ind)));
MY = MX;
MX = [MX, XR(bound_red(ind))];
MY = [MY, YR(bound_red(ind))];
% MX = [XO(bound_orange(ind)).^2, YO(bound_orange(ind)).^2, XO(bound_orange(ind)), YO(bound_orange(ind)), XO(bound_orange(ind)).*YO(bound_orange(ind)), XO(bound_orange(ind))*0+1, XR(bound_red(ind))];
% MY = [XO(bound_orange(ind)).^2, YO(bound_orange(ind)).^2, XO(bound_orange(ind)), YO(bound_orange(ind)), XO(bound_orange(ind)).*YO(bound_orange(ind)), XO(bound_orange(ind))*0+1, YR(bound_red(ind))];
MX = rref(MX);
MY = rref(MY);


o2rx = MX(:,end);
o2ry = MY(:,end);
plot(XR(bound_red),YR(bound_red));
plot(XO(bound_orange),YO(bound_orange));

% Idea that averaging all combination of several points at a 6 feature
% vector shot will give the best coverage parameters for all data


% average the result
% o2rx = mean(o2rx,2);
% o2ry = mean(o2ry,2);
% vector = [XO.^2, YO.^2, XO, YO,XO.*YO, XO*0+1];
vector = xy_feature(XO,YO);
XOP = o2rx.'*vector.';
YOP = o2ry.'*vector.';

plot(XOP,YOP,'.g')



figure
ind = xf < split;
xfr = xf(ind);
yfr = yf(ind);
xfo = xf(logical(1-ind));
yfo = yf(logical(1-ind));

% vector = [xfo.^2, yfo.^2, xfo, yfo,xfo.*yfo, xfo*0+1];
vector = xy_feature(xfo,yfo);
xo2r = o2rx.'*vector.';
yo2r = o2ry.'*vector.';
plot(xfr,yfr,'.b')
hold on
plot(xo2r,yo2r,'.r')
hold off
axis equal

% average points based on several frames of the same latex bead
% added = (1:numel(XO))*0;

% Rearrange position arrays so indices match for bead images in respective
% channels
for i = 1:numel(XO)
    dist = ((XOP-XR(i)).^2 + (YOP-YR(i)).^2 ).^0.5;
    ind = find(dist == min(dist));
    indexs(i,1) = ind;
end
XOP =XOP(indexs);
YOP =YOP(indexs);
% ind = [3, 8,4, 13];
% % ind = 4; 11
% XT = XT(:);
% YT = YT(:);
% XRT =XRT(:);
% YRT =YRT(:);
% XT = [XT(1:end-4);XOP(ind).'];
% YT = [YT(1:end-4);YOP(ind).'];
% XRT= [XRT(1:end-4);XR(ind)];
% YT = [YRT(1:end-4);XR(ind)];
% % XO = x;
% % YO = y;
% close all
% MX = xy_feature(XT,YT);
% MY = MX;
% MX = [MX, XRT];
% MY = [MY, YRT];
% MX = [XO(bound_orange(ind)).^2, YO(bound_orange(ind)).^2, XO(bound_orange(ind)), YO(bound_orange(ind)), XO(bound_orange(ind)).*YO(bound_orange(ind)), XO(bound_orange(ind))*0+1, XR(bound_red(ind))];
% MY = [XO(bound_orange(ind)).^2, YO(bound_orange(ind)).^2, XO(bound_orange(ind)), YO(bound_orange(ind)), XO(bound_orange(ind)).*YO(bound_orange(ind)), XO(bound_orange(ind))*0+1, YR(bound_red(ind))];
% MX = rref(MX);
% MY = rref(MY);


% o2rx = MX(:,end);
% o2ry = MY(:,end);


% Data representation
figure
plot(xf,yf,'.k')
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Composite bead data');
saveas(gcf,'Scatter_data.png');


hold off
plot(xfr,yfr,'.r');
hold on
plot(xfo,yfo,'.b');
plot([split,split],[min(yfr),max(yfr)],'k')
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Composite bead data');
saveas(gcf,'Scatter_data_with split.png');
% hold off


plot(XR(bound_red),YR(bound_red),'r','LineWidth', 2);
hold on
plot(XO(bound_orange),YO(bound_orange),'b','LineWidth', 2);
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Composite bead data');
saveas(gcf,'Scatter_data_with boundary.png');
hold off

xo2r = o2rx.'*vector.';
yo2r = o2ry.'*vector.';
plot(xo2r,yo2r,'.b')
hold on
plot(xfr,yfr,'.r')
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Composite bead data');
legend('Orange Beads','Red Beads')
saveas(gcf,'Corrected_data boundary.png');

xlim([54, 65])
ylim([119, 128])
saveas(gcf,'Corrected_data zoom.png');
% 
xlim([60.4, 61.1])
ylim([124.1, 124.6])
saveas(gcf,'Corrected_data bad_spot super zoom.png');

xlim([41.15, 41.65])
ylim([161.6, 162])
saveas(gcf,'Corrected_data_good zoom.png');
hold off
[xf,yf] = meshgrid(min(xfo):(max(xfo)-min(xfo))/20:max(xfo),min(yfo):(max(yfo)-min(yfo))/20:max(yfo));
xfo = xf(:);
yfo = yf(:);
% vector = [xfo.^2, yfo.^2, xfo, yfo,xfo.*yfo, xfo*0+1];
vector = xy_feature(xfo,yfo);
xo = o2rx.'*vector.';
yo = o2ry.'*vector.';
% plot(xfr,yfr,'.r');
% hold on
figure
plot(xfo,yfo,'.r')
hold on
plot(xo,yo,'.b');
legend('Theoretical grid','Resulting transform')

% xlim([30, 160])
% ylim([60, 170])
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Theoretical grid');
saveas(gcf,'regular_sized grid.png');
save('2_color_calibration.mat','o2rx','o2ry', 'split');
figure
plot(xo2r,yo2r,'.b')
hold on
plot(xfr,yfr,'.r')
axis equal
xlabel('X position in pixels')
ylabel('Y position in pixels')
title('Composite bead data');
legend('Orange Beads','Red Beads')