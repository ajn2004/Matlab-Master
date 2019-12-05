function i1 = get_psf(x0,y0,sx,sy)
pixw = 9;
[xpix,ypix] = meshgrid(-pixw:pixw);
i1x = 1/2.*(erf((xpix - x0 + 1/2)./(2*sx^2)^0.5)-erf((xpix - x0 - 1/2)./(2*sx^2)^0.5)); % error function of x integration over psf for each pixel
i1y = 1/2.*(erf((ypix - y0 + 1/2)./(2*sy^2)^0.5)-erf((ypix - y0 - 1/2)./(2*sy^2)^0.5)); % error function of y integration over psf for each pixel
i1 = i1x.*i1y;
i1 = i1./sum(i1(:));