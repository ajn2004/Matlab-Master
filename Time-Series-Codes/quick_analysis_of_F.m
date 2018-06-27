%% Quick analysis of F
% This script will allow the user to select a number of points and create a
% series of images showing the average F over a user selected number of
% boutons as well as a where cutoffs are to measure the noise associated
% with that signal

clearvars
close all
f = figure('Name','F Trace Analysis');
try
    load('back_subtract.mat');
catch mi1 = 0;
end
% c = scrub_config();
% i1 = (readtiff()-mi1)/em_gain(c.Gain);
i1 = (readtiff()-mi1)/em_gain(300);
[m,n,o] = size(i1);
pixw = 7;
tex = 0.0309;
stim = 100;
stims = 1;
str = 10;
si1 = std(i1,1,3);
imagesc(si1)
exps = 141;
num = input('How many points?');
[x,y] = ginput(num);

wind = -pixw:pixw;

for i = 1:num
    sub1 = i1(round(y) + wind, round(x) + wind,:);
    mfl = sum(sum(sub1));
    sfluor = sum(sum(sub1));
    sfluor = sfluor(:);
    mfluor(:,i) = mfl(:);
end
% mfluor = mean(mfluor,2);
% bkg = mean(mfluor(68:stim-1));
% nois = std(mfluor(68:stim-1));
t = (1:o)*tex;
% df = mfluor - bkg;
% ndf = df/bkg;
% nf = mfluor/bkg;
tbgrp = uitabgroup(f);
fl = uitab(tbgrp,'Title','F and stims');
tg = uitabgroup(fl);
ax1 = axes(uitab(tg,'Title','1'));
plot(ax1,t,mfluor);
hold on
ffluor = mfluor(exps:end);
for i = 1:stims
    plot([stim*tex + (i-1)/str,stim*tex + (i-1)/str],[min(mfluor),max(mfluor)],'r');
end
xlabel(ax1,'Seconds')
ylabel(ax1,'F (AU)');
title('HILO vGlut 10 Stim 10Hz')
% title(ax1,'non-HILO vGlut 10 Stim 10Hz')
legend(ax1,'F','Stims')


%
% ax2 = axes(uitab(tg,'Title','2'));
% % close all
% plot(ax2,t,mfluor)
% sf = gausssmooth(mfluor,20,50);
% hold on
% plot(ax2,t,sf);
% legend(ax2,'F','Smoothed F');
% % title(ax2,'non-HILO vGlut 10 Stim 10Hz')
% title('HILO vGlut 10 Stim 10Hz')
% xlabel(ax2,'Seconds')
% ylabel(ax2,'F (AU)');
%
%
% ax = axes(uitab(tg,'Title','3'));
% plot(ax,t,mfluor)
% sf = gausssmooth(mfluor,20,50);
% hold on
% plot(t,sf);
% plot([5.778, 5.778],[min(mfluor), max(mfluor)],'r');
% plot([2.101, 2.101],[min(mfluor), max(mfluor)],'r');
% plot([stim*tex, stim*tex],[min(mfluor), max(mfluor)],'r');
% legend('F','Smoothed F','Cutoff');
% % title('non-HILO vGlut 10 Stim 10Hz')
% title('HILO vGlut 10 Stim 10Hz')
% xlabel('Seconds')
% ylabel('F (AU)');
% ind = t > 5.778;
% mf = mfluor(ind) - sf(ind);
%
%
% ax = axes(uitab(tg,'Title','4'));
% plot(t(ind),mf);
% title('Ave Subtracted post-stim');
% % title('HILO vGlut 10 Stim 10Hz')
% xlabel('Seconds')
% ylabel('F (AU)');
%
% ax = axes(uitab(tg,'Title','dF/F'));
% plot(t,(mfluor-bkg)/bkg);
% title('Ave Subtracted post-stim');
% % title('HILO vGlut 10 Stim 10Hz')
% xlabel('Seconds')
% ylabel('F (AU)');
% snr = (max(mfluor)-bkg)/std(mf)
% % close all
%

