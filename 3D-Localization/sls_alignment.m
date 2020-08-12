% function diff = sls_alignment(sig1, sig2)
% This function is to return the index of most overlap between two sigma
% curves by finding the frame intercept in the frame v sigma curve then
% comparing those intercepts.
% Each sig input should be packed in the form of [fnum, sx,sy]
% AJN 
% Unpacking
% Smooth inputs to avoid noise w/ spline fit
sig2 = sigc;
sxy1 = gausssmooth(sig1(:,2).^2 - sig1(:,3).^2,4,10);
sxy2 = gausssmooth(sig2(:,2).^2 - sig2(:,3).^2,4,10);
f1 = sig1(:,1);
f2 = sig2(:,1);
% fout = min([f1,f2]):max([f1,f2]);
close all

plot(f1, sxy1)
hold on
plot(f2,sxy2)
% 
hold off
% figure
% [c, lags] = xcorr(sxy2,sxy1);
cost = [];
displacements = [ -80:80];
for j = displacements
    cost = [cost; scroll_cost(f1,sxy1,f2-j,sxy2)];
end
plot(displacements,cost)
figure
offset = displacements(cost == min(cost));
% stem(lags,c)
% offset = lags(find(c == max(c)));
% figure
plot(f1,sxy1);
hold on
title('corrected')
plot(f2 - offset, sxy2);
hold off

function cost = scroll_cost(x1,y1,x2,y2)
% Find our 'bounds' 
% For this problem we're dealing with frames which are integers that correspond
% to our index of a sort
lower_bound = max([min(x1),min(x2)]); % we want to take the greater of the two min values, this ensures overlap
upper_bound = min([max(x1),max(x2)]); % Taking the minimal upper bound ensures overlap
cost = 0;
for k = lower_bound:upper_bound
    ind_1 = find(x1 == k);
    ind_2 = find(x2 == k);
    if ~isempty(ind_1) && ~isempty(ind_2)
        cost = cost + (y1(ind_1) - y2(ind_2))^2;
    end
end
end



