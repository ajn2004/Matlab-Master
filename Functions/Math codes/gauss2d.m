function z = gauss2d(vec,pixw)
[xpix,ypix] = meshgrid(-pixw:pixw); % Create mesh grids for desired window
i1 = xpix.*0;
x0 = vec(1);
y0 = vec(2);
N = vec(3);
sx = vec(4);
sy = vec(5);
B = vec(6);
% Create a gaussian
i1x = 1/2.*(erf((xpix - x0 + 1/2)./(2*sx^2)^0.5)-erf((xpix - x0 - 1/2)./(2*sx^2)^0.5)); % error function of x integration over psf for each pixel
i1y = 1/2.*(erf((ypix - y0 + 1/2)./(2*sy^2)^0.5)-erf((ypix - y0 - 1/2)./(2*sy^2)^0.5)); % error function of y integration over psf for each pixel
z = N * i1x.*i1y+B;