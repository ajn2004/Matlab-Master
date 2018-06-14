function [i2] = chop_up(i1,pixw,fnum, cents, aves)
% chop up an image given frame number and central pixel
[m,n,o] = size(i1);
i2 = [];
wind = -pixw:pixw;
for j = 1:numel(fnum) % loop over all frames
    
    i2 = cat(3,i2,i1(cents(j,1) + wind, cents(j,1) + wind,fnum(j)-aves:fnum(j)-1));

end

