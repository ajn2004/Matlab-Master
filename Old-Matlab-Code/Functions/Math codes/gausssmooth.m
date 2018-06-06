function s3 = gausssmooth(s1, gwidth, wind)

%% Gaussian for smoothing
try
s2 = [ones(1,wind)*s1(1), s1, s1(end)*ones(1,wind)];
catch
    s2 = [ones(wind,1)*s1(1); s1; s1(end)*ones(wind,1)];
end
x = -wind:wind;
gauss = exp(-x.*x/(2*gwidth^2));
ngauss = gauss./sum(gauss);
ssmooth = conv(s2(:),ngauss(:),'same'); % smooth g
s3 = ssmooth(wind+1:end-wind);
try
    polyfit(s1,s3,1);
catch
    s3=s3.';
end
end