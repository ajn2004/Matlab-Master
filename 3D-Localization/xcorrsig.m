function diff = xcorrsig(sig1, sig2)
% This function is to return the index of most overlap between two sigma
% curves by finding the frame intercept in the frame v sigma curve then
% comparing those intercepts.
% Each sig input should be packed in the form of [fnum, sx,sy]
% AJN 
% Unpacking
% Smooth inputs to avoid noise w/ spline fit
sxy1 = gausssmooth(sig1(:,2).^2 - sig1(:,3).^2,4,10);
sxy2 = gausssmooth(sig2(:,2).^2 - sig2(:,3).^2,4,10);
f1 = sig1(:,1);
f2 = sig2(:,1);
% fout = min([f1,f2]):max([f1,f2]);

% Spline interpolation onto the output

wind = -5:5;
mid1 = round(numel(sxy1)/2);
mid2 = round(numel(sxy2)/2);

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



