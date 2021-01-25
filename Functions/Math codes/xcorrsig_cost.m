function cost = xcorrsig_cost(x1,y1,x2,y2)
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