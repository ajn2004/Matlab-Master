%%% Building an analysis for localizing pHluorin molecules
clearvars
close all
clc
%% User Variabls
prescans = 100; % number of 'prescanned' frames
stims = 98; % stimulus interval
tex = 1/32.625;% total frame period
pixw = 7;
%% Non-user Section
load('back_subtract.mat'); % load background
i1 = (readtiff() - mi1)/33.33; % puts the image into photons without dark current
imagesc(std(i1,1,3))
[x,y] = ginput(1);
wind = -pixw:pixw;
hold on
plot([x - pixw,x + pixw, x+pixw, x-pixw,x-pixw],[y-pixw, y-pixw, y+pixw, y+ pixw,y-pixw],'r');
hold off
sub = i1(round(y) +wind, round(x) + wind,:);
mf = mean(mean(sub));
mf = mf(:);
figure
t = (1:numel(mf))*tex;
plot(t,mf);
hold on
% Calculating Stims

