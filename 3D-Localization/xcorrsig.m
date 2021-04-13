function diff = xcorrsig(sig1, sig2)
% This function is to return the index of most overlap between two sigma
% curves by finding the frame intercept in the frame v sigma curve then
% comparing those intercepts.
% Each sig input should be packed in the form of [fnum, sx,sy]
% AJN 
% Unpacking
% Smooth inputs to avoid noise w/ spline fit

% sxy1 = gausssmooth(sig1(:,2).^2 - sig1(:,3).^2,4,10);
% sxy2 = gausssmooth(sig2(:,2).^2 - sig2(:,3).^2,4,10);
f1 = sig1(:,1);
f2 = sig2(:,1);
n = numel(f1);
sx1 = sig1(:,2);
sy1 = sig1(:,3);
sx2 = sig2(:,2);
sy2 = sig2(:,3);
frame1 = sig1(:,1);
frame2 = sig2(:,1);
% [c, lags] = xcorr(sxy2,sxy1);
cost = [];
displacements = [-floor(n/4):floor(n/4)];
for j = 1:numel(displacements)
    dist = 0;
    for i = 1:numel(sx2)
        ind = frame1 == frame2(i) + displacements(j);
        if sum(ind) == 1
            dist = dist + (sx1(ind)-sx2(i)).^2 + (sy1(ind)-sy2(i))^2;
        end
    end
    cost = [cost; dist];
%     cost = [cost; xcorrsig_cost(f1,sxy1,f2-j,sxy2)];
end
s_cost = gausssmooth(cost,4,10);
diff = displacements(s_cost == min(s_cost(s_cost > 0)));
diff = diff(1);
% stem(lags,c)
% offset = lags(find(c == max(c)));
% figure


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


end
