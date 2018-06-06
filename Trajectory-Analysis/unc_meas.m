%% Uncertainty Measurer
clearvars
close all
clc;

% [fname, fpath] = uigetfile('*traj.mat');
% load([fpath,fname]);
files = dir('*traj.mat');

dx = [];
dy = [];
dz = [];
for j = 1:numel(files)
    load(files(j).name,'xf_all','yf_all','zf_all', 'trajec','q');
%     zf_all = getdz(sigx_all,sigy_all)/q;
trs = numel(trajec);
% for i = 1:trs
%     dx = [dx;1000*q*(xf_all(trajec(i).t) - mean(xf_all(trajec(i).t)))];
%     dy = [dy;1000*q*(yf_all(trajec(i).t) - mean(yf_all(trajec(i).t)))];
%     dz = [dz;1000*q*(zf_all(trajec(i).t) - mean(zf_all(trajec(i).t)))];
% end
trs = numel(trajec);
for i = 1:trs
    dx = [dx;1000*q*diff(xf_all(trajec(i).t))];
    dy = [dy;1000*q*diff(yf_all(trajec(i).t))];
    dz = [dz;1000*q*diff(zf_all(trajec(i).t))];
end
end
figure
hx = histogram(dx);
lx = (sum(dx.^2)/(numel(dz)-1))^0.5
ly = (sum(dy.^2)/(numel(dz)-1))^0.5
lz = (sum(dz.^2)/(numel(dz)-1))^0.5
hold on
hy = histogram(dy,'Normalization','Probability');
legend('X','Y');
xlabel('Distance from mean position in nm');
ylabel('Probability');
title('Lateral Uncertainty');
hold off
figure
hold on
histogram(dz,'Normalization','Probability');
histogram(dx,'Normalization','Probability');
histogram(dy,'Normalization','Probability');
legend('Z','X','Y');
ylabel('Probability');
xlabel('Distance from mean position in nm');
ylabel('Probability');
title('Axial Uncertainty');
hold off
% legend('dx','dy','dz');