function [Cnew, clus_new] = changerows(Cold,j, clus_id)
% this is intended to be a recrusive function
Cnew = Cold;
clus_new = clus_id;

for ii = 1 :numel(Cnew(1,:))
    if Cnew(j,ii) == 1 && clus_new(ii) == 0
        clus_new(ii) = clus_new(j);
        [Cnew, clus_new] = changerows(Cnew, ii, clus_new);
    end
end