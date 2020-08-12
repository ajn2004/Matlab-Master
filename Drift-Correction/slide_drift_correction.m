% Trial Drift Correction Idea
% Attempt to minimize variation of a data set by iteratively moving
% partitions of localization data in x,y and z respectively

clearvars;
close all;
clc;

fms = 1000;
fpath = 'C:\Users\AJN Lab\Dropbox\Data\6-3-19 cono-halo Neurons\Analysis\toleranced\DC\';
fname = 'Cell4_dz10_r1_dast_tol_dc.mat';
load([fpath,fname],'ncoords','xf_fixed','yf_fixed','q','framenumber');
r = str2num(fname(strfind(fname,'_r')+2));
zf = func_shift_correct(ncoords(:,3)*q,framenumber,r);
xf = ncoords(:,1)*q;
yf = ncoords(:,2)*q;
scanr = -1:0.005:1;
chunks = floor((framenumber(end) - 1)/fms) + 1;
ind0 = framenumber <= fms+1; % initialize a data set
% plot3(xf(ind0),yf(ind0),zf(ind0),'.')
hold on
for i = 2:chunks
    fm1 = (i-1)*fms + 1;
    fm2 = i*fms + 1;
    ind1 = framenumber >= fm1 & framenumber <=fm2;
%     plot3(xf(ind1),yf(ind1),zf(ind1),'.')
    for j = 1:numel(scanr)
        testx = [xf(ind0);xf(ind1) - scanr(j)];
        testy = [yf(ind0);yf(ind1) - scanr(j)];
        testz = [zf(ind0);zf(ind1) - scanr(j)];
        stdx(j) = std(testx);
        stdy(j) = std(testy);
        stdz(j) = std(testz);
%         stdx(j) = mean(xf(ind0)) - (mean(xf(ind1)) - scanr(j));
%         stdy(j) = mean(yf(ind0)) - (mean(yf(ind1)) - scanr(j));
%         stdz(j) = mean(zf(ind0)) - (mean(zf(ind1)) - scanr(j));
    end
    indx = find(stdx == min(stdx));
    indy = find(stdy == min(stdy));
    indz = find(stdz == min(stdz));
    driftx(i) = scanr(indx);
    drifty(i) = scanr(indy);
    driftz(i) = scanr(indz);
%     plot(scanr,stdx)
%     hold on
    plot3(xf(ind1)-sum(driftx(2:i)),yf(ind1)-sum(drifty(2:i)),zf(ind1)-sum(driftz(2:i)),'.')
    hold on
    ind0 = ind1;
end
axis equal
hold off
% figure
% plot(drift)