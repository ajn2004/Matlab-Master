function [dps] = cpu_peaks(i1, thrsh,pixw)

[m,n,o] = size(i1); % get size of input image
dps = zeros(m,n,o); % preallocate output array
wind = -pixw:pixw; % define a window to find the highest point around
tsh = i1.*(i1>=thrsh);

for i = 1:o % loop over all frames
%     tsh(:,:,i) = i1(:,:,i).*(i1(:,:,i)>=(max(max(i1(:,:,i)))*thrsh/100));
    [row,col] = find(tsh(:,:,i) > 0); % find all pixels above threshold
    for j = 1:numel(row) % loop over all pixels above threshold in frame
        if row(j)+wind(1) > 0 && row(j)+wind(end) < m && col(j)+wind(1) > 0 && col(j)+wind(end) < n
        sub = i1(row(j)+wind,col(j)+wind,i); % grab subbed region around the threshold point
        s2t = sub.*(sub>=thrsh*0.5);
        if sum(s2t(:)) >=2 && tsh(row(j),col(j),i) == max(sub(:)) % if current pixel is the maximal pixel in the region
            dps(row(j),col(j),i) = 1; % if condition is satisfied assign a 1 to the dps variable location
        end
        end
    end
end
    