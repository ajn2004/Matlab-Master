% compare bkgn subtraction algorithms
clearvars -except i1
close all
if ~exist('i1')
    i1 = readtiff();
end

i1 = i1(:,:,1);

irb = rollingball(i1);
ibp = bandpass(i1);
figure('units','Normalized','Outerposition',[0 0 1 1])
subplot(3,3,1);
imagesc(i1);
title('Raw')
axis image
subplot(3,3,2);
C = xcorr2(i1,i1);
imagesc(C/max(C(:)));
title('xcorr')
axis image
subplot(3,3,3);
imagesc(i1-i1);
title('"diff"')
axis image
subplot(3,3,4);
imagesc(irb);
title('RollingBall')
axis image
subplot(3,3,5);
Crb = xcorr2(i1,irb);
imagesc(Crb/max(Crb(:)));
axis image
subplot(3,3,6);
imagesc(i1-irb);
axis image
subplot(3,3,7);
imagesc(ibp);
title('Bandpass')
axis image
subplot(3,3,8);
Cbp = xcorr2(i1,ibp);
imagesc(Cbp/max(Cbp(:)));
axis image
subplot(3,3,9);
imagesc(i1-ibp);
axis image