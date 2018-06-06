%% Trajectory Spinner
% This script will create a figure from localized data to spin trajectory
% information around a defined theta and phi path then record the resulting
% frames
[fname, fpath] = uigetfile('*traj.mat');
load([fpath,fname]);
Trajec_viewer();
view([0,90])
[x,y] = ginput(2);
xlim([x(1),x(2)])
ylim([y(2),y(1)])
zlim([-1,1])
title('Rotated Bouton');
xlabel('Position [um]');
ylabel('Position [um]');
zlabel('Position [um]');
view([0,90]);
M(1) = getframe(gcf);
eles=(90:-1:0);

for i = 1:numel(eles)
    view([0,eles(i)]);
    M(numel(M)+1) = getframe(gcf);
end

azs = 1:360;
for i = 1:numel(azs)
    view([azs(i),0]);
    M(numel(M)+1) = getframe(gcf);
end