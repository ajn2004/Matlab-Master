%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim viewer
%
% View the results of a sim
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;
load('recent_sim.mat');

% for i = 1:op
    xd = (xfa - xf0);
    sxd = sum(xd.^2,1)./(numel(xd(:,1)) - 1);
    mx = mean(sxd,1);
    yd = (yfa - yf0);
    syd = sum(yd.^2,1)./(numel(yd(:,1)) - 1);
    my = mean(syd,1);
    zd = (zfa - zf0);
    szd = sum(zd.^2,1)./(numel(zd(:,1)) - 1);
    mz = mean(szd,1);
% end

% Build gaussians

plot(Ns,mx*128);
hold on
plot(Ns,my*128);
hold off
legend('X-unc','Y-unc');
set(gca,'Xscale','log');
set(gca,'Yscale','log');
xlabel('log(photons)');
ylabel('log(nm)');
title('Uncertainty in X and Y vs fitted N')
figure
plot(Ns,mz*1000);
set(gca,'Xscale','log');
set(gca,'Yscale','log');

a = polyfit(log(Ns(5:end)),log(mz(5:end)*1000),1);
hold on
plot([150 150],[5 40],'r');
plot([350 350],[5 40],'r');
plot([1000 1000],[0.1 10],'g');
plot([4000 4000],[0.1 10],'g');
ylim([0.09 1000]);
xlim([100 7000]);
xlabel('log(photons)');
ylabel('log(nm)');
title('Axial Uncertainty vs Photons Measured');
hold off