% Release probability Estimator
% This is a code to estimate the release probability of 'local' data set.
% It will do so initially by measuring a region's fluoresence in the first
% background frames and comparing that to the post-stimulus frames
close all; clearvars; clc;

%% Load data into variable A
files = dir('*tif');
A = {};
f = figure('Units','Normalized','Outerposition',[0 0 1 1])
h = uitabgroup(f);
for i = 1:numel(files) % This section builds an array of idecies to keep track of chronology
    try
        ind(i) = str2num(files(i).name(7:end-4));
    catch
        ind(i) = 0;
    end
end
[B,I] = sort(ind); % I gives least to greatest order for ind variable
i1 = [];
for i = 1:numel(files)
    A{i} = readtiff(files(I(i)).name); % load in order of I to a cell variable
    i1 = cat(3,i1,A{i}); % Build an image variable for quick identification
end

% Image selection and segmentation
him = uitab(h,'Title','Max Projection');
ax = axes(him);
imagesc(ax,max(i1,[],3));
axis image
colormap('jet')
[x,y] = ginput(1);
x = round(x);
y = round(y);
pixw = 5;
wind = -pixw:pixw;
hold on
plot(ax,x+wind,(y+wind)*0+y+pixw,'m')
plot(ax,x+wind,(y+wind)*0+y-pixw,'m')
plot(ax,(x+wind)*0+x+pixw,y+wind,'m')
plot(ax,(x+wind)*0+x-pixw,y+wind,'m')
% Quick description of the region
subi = i1(y + wind, x + wind,:);
[m,n,o] = size(subi);
sfluor  = reshape(sum(sum(subi)),o,1);
mf = uitab(h,'Title','F-Trace');
ax = axes(mf);
plot(ax,sfluor)
xlabel(ax,'Frame')
ylabel(ax,'Photons summed over the region');

% add stims in the region
fps = 6; % frames per set
stim = 4; % frame stim occurs on in the set
stim = stim - 0.5; % offset stim to halfway between frames

hold on
for i = 1:numel(A)-1 % add stims
    plot(ax,[stim+i*fps, stim+i*fps],[min(sfluor), max(sfluor)],'r')
    plot(ax,[(i-1)*fps+0.5,0.5+(i-1)*fps],[min(sfluor), max(sfluor)],'k')
end
legend(ax,'Signal','Stimulus','Separation')
hold off
% More detailed look at individual stimulus
% figure
hold on
i3 = zeros(2*pixw+1);
for i = 1:numel(A)
    i2 = A{i}(y + wind, x + wind,:);
    sf2(i,:) = reshape(sum(sum(i2)),1,fps);
    sf2(i,:) = sf2(i,:) - mean(sf2(i,:));
    
    df2(i) = mean(sf2(i,4:6))-mean(sf2(i,1:3));
    if df2(i) > 0
        i3 = i3 + (i2(:,:,4) - mean(i2(:,:,1:3),3));
    end
end
i3 = i3/sum(df2>0);
hd = uitab(h,'Title','Hist Df2');
ax = axes(hd);
histogram(ax,df2,'BinWidth',10)
thresh  = input('What should the threshold be?')
ind = df2 > thresh;
bt = uitab(h,'Title','Box Plots');
bg = uitabgroup(bt);
rt = uitab(bg,'Title','Releases');
nrt = uitab(bg,'Title','Non-Releases');
ax = axes(rt);
boxplot(ax,sf2(ind,:),1:6)
title(ax,'Probable release')
ax = axes(nrt)
boxplot(ax,sf2(logical(1-ind),:),1:6)
title(ax,'Probable non-release')
[h1,p] = ttest2(sf2(ind,1),sf2(logical(1-ind),1))
[h2,p] = ttest2(sf2(ind,2),sf2(logical(1-ind),2))
[h3,p] = ttest2(sf2(ind,3),sf2(logical(1-ind),3))
[h4,p] = ttest2(sf2(ind,4),sf2(logical(1-ind),4))
[h5,p] = ttest2(sf2(ind,5),sf2(logical(1-ind),5))
[h6,p] = ttest2(sf2(ind,6),sf2(logical(1-ind),6))
sum(df2>thresh)/numel(df2)
ext = uitab(h,'Title','Example Average Difference');
ax = axes(ext)
surf(ax,i3);
axis image
 