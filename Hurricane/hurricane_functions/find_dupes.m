function ind = find_dupes(cents,fnum)
% function for finding duplicates and returning a variable indexing all
% redundant copies of that variable
ind = [];
for i = 1:numel(fnum) % loop through all found molecules
    if ~ismember(i,ind) % if current molecule isn't already a duplicate continue
        indy = []; % clear temp index variable
        indy = find(fnum == fnum(i) & cents(:,1) == cents(i,1) & cents(:,2) == cents(i,2));
        % find all indices whose paramters match those of current ith index
        if numel(indy) > 1 % if there are more than just 1 copy
            ind = [ind;indy(2:end)]; % add list of copies to overall index
        end
    end
end