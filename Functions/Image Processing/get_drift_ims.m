function [drift] = get_drift_ims(ims)
% a function to get the drift between images using the cross correlation
% method

[m,n,o] = size(ims); % find number of images in data set
wind = -5:5; % look at 9x9 window centered around the maximum correlation value
[X,Y] = meshgrid(wind);
drift = [0,0]; % no drift on the first frame
irb = rollingball(ims);
ibp = bandpass(ims);
figure('Units','Normalized','OuterPosition',[0 0 1 1]);
for i = 2:o % run over every frame in data set
    C = xcorr2(ims(:,:,i-1),ims(:,:,i)); % x correlation
    Crb = xcorr2(irb(:,:,i-1),irb(:,:,i)); % x correlation
    Cbp = xcorr2(ibp(:,:,i-1),ibp(:,:,i)); % x correlation
    subplot(1,3,1);
    imagesc(C)
%     % perform calculation for raw data, rolling ball subtracted, and band
%     % passed data
    [row,col] = find(C == max(C(:))); % find peak
    hold on
    plot(col,row,'rx');
    hold off
    axis image
%     drawnow
%     
%     subC = C(row+wind, col + wind); % segment area around peak
%     subC = subC/max(subC(:)); % normalize peak value to 1
%     [ang] = get_elip_ang(subC,3,4);
%     [fits] = func_mle_crlb(subC,0,0,3,ang); % perform a gaussian fit to image
%     clear C subC row col
%     Crb = rollingball(Crb);
    [rrb,crb] = find(Crb == max(Crb(:))); % find peak
    subC = Crb(row+wind, col + wind); % segment area around peak
    xcm = sum(sum(X.*subC))/sum(subC(:));
    ycm = sum(sum(Y.*subC))/sum(subC(:));
    rrb = rrb + ycm;
    crb = crb + xcm;
%     subC = subC/max(subC(:)); % normalize peak value to 1
%     [ang] = get_elip_ang(subC,3,4);
%     [frb] = func_mle_crlb(subC,0,0,3,ang); % perform a gaussian fit to image
%     clear C subC row col
subplot(1,3,2);
imagesc(Crb);
hold on
plot(crb,rrb,'rx');
axis image
hold off
%     Cbp = rollingball(Cbp);
    [rbp,cbp] = find(Cbp == max(Cbp(:))); % find peak
    
    subC = Cbp(row+wind, col + wind); % segment area around peak
    xcm = sum(sum(X.*subC))/sum(subC(:));
    ycm = sum(sum(Y.*subC))/sum(subC(:));
    rbp = rbp + ycm;
    cbp = cbp + xcm;
%     subC = subC/max(subC(:)); % normalize peak value to 1
%     [ang] = get_elip_ang(subC,3,4);
%     [fbp] = func_mle_crlb(subC,0,0,3,ang); % perform a gaussian fit to image
    subplot(1,3,3);
    imagesc(Cbp)
    hold on
    plot(cbp,rbp,'rx');
    axis image
    hold off
    drawnow
    % average the final results
%     drift = [drift;drift(i-1,1) + col + (fits(1) + fbp(1))/2-n,drift(i-1,2) + row + (fbp(2) + fits(2))/-m];
    drift = [drift;drift(i-1,1) + (crb(1)+cbp(1))/2 - n, drift(i-1,2) + (rrb(1)+rbp(1))/2-m];
%     drift = [drift;(crb(1)+cbp(1))/2 - n, (rrb(1)+rbp(1))/2-m];
%     drift = [drift; (fits(1)+frb(1)+fbp(1))/3, (fits(2)+frb(2)+fbp(2))/3]; % record value in output variable
    clear C subC fits fbp frb crlbs llv row col % clear values that are used in the calculation
    M(i-1) = getframe(gcf);
end
movie2gif(M,'xcorrelation.gif','DelayTime',0.05,'LoopCount',Inf);