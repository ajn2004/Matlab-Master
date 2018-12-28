function diff = xcorrsig2(sig1, sig2)
% This function is to return the index of most overlap between two sigma
% curves by finding the frame intercept in the frame v sigma curve then
% comparing those intercepts.
% Each sig input should be packed in the form of [fnum, sx,sy]
% AJN 
% Unpacking
% Smooth inputs to avoid noise w/ spline fit
sx1 = gausssmooth(sig1(:,2),4,10);
sy1 = gausssmooth(sig1(:,3),4,10);
sx2 = gausssmooth(sig2(:,2),4,10);
sy2 = gausssmooth(sig2(:,3),4,10);


f1 = sig1(:,1);
f2 = sig2(:,1);
% fout = min([f1,f2]):max([f1,f2]);

% Spline interpolation onto the output

wind = -15:15;
mid1 = round(numel(sxy1)/2);
mid2 = round(numel(sxy2)/2);

for i = 1:(numel(sx2)-1)
        inds = i:(numel(sx2)-1);
        cost(i) = sum((sx1(inds) - sx2(inds)).^2 + (sy1(inds) - sy2(inds)).^2)/numel(inds);
        plot(sx1(inds))
        hold on
        plot(sx2(inds))
        drawnow
end
ssx1 = sxy1(mid1+wind);
try
ssx2 = sxy2(mid2+wind);
fs1 = f1(mid1+wind);
fs2 = f2(mid2+wind);

a = polyfit(ssx1,fs1,1);
b = polyfit(ssx2,fs2,1);

diff = round(a(2)-b(2));
catch lsterr
    diff = 0;
end



