% Comparing filters
[m,n,o] = size(dip1);
% M = [];
clear M
rip1 = rollingball(dip1);
mins = (0.61*0.51/1.3)/(2*q);
maxs = 2*mins;
dip3 = bandpass(dip1,mins,maxs);
lip2 = lp_filt(dip1,30);
rlims = [min(min(min(rip1))), max(max(max(rip1(:,:,1:2))))];
blims = [min(min(min(dip3))), max(max(max(dip3(:,:,1:2))))];
clims = [min(min(min(dip1))), max(max(max(dip1(:,:,1:2))))];
flims = [min(min(min(lip2))), max(max(max(lip2(:,:,1:2))))];
for i = 1:o
    subplot(2,2,1);
    imagesc(dip1(:,:,i),clims);
    axis image
    title('Raw Difference Frame')
    subplot(2,2,2);
    imagesc(rip1(:,:,i),rlims);
    axis image
    title('Rollingball Subtraction of Difference Image')
    subplot(2,2,3);
    imagesc(dip3(:,:,i),blims);
    axis image
    title('Bandpass filter');
    subplot(2,2,4);
    imagesc(lip2(:,:,i),flims);
    title('Low Pass filter');
    axis image
    pause(0.1)
        M(i) = getframe(gcf);

end