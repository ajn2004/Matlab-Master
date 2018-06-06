function i2 = ave_subtract(i1,stim)
% This function will average the stim previous frames and subtract them
% from the stim frame
for i = 1:floor(numel(i1(1,1,:))/stim)
    i1sub = mean(i1(:,:,(i-1)*stim +1:i*stim-1),3);
    i2(:,:,i) = i1(:,:,i) - i1sub;
end
end