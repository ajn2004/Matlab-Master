function [x1, y1 ] = func_build_gauss(x,y,nmin, nmax, sigmin, sigmax)
% Builds a series of gaussian images based on inputs x and y for neural
% learning
x1= x;
y1= y;
m = numel(y);
n = numel(x(1,:))^0.5;

w = (n-1)/2;
rbox = w;
[xpix, ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
% xpix = gpuArray(xpix);
% ypix = gpuArray(ypix);
wbox = 2*rbox+1;
i1 = xpix.*0;

% variables of image
B = 0.5;
% ntrue = 1000;
x0true = 0;
y0true = 0;
sigma2 = 1.5;   %this would be rayleigh radius in pixel space
sigma = sigma2/2;   % sigma used for fitting
% r0t_pix = 2*sigma; % 1/e^2 radius
r0t_um = .2;        % Rayleigh radius in um
% pix2pho = 1;
q = r0t_um / (sigma2);
% wbox_um = q*wbox;

% frames = 1000;
% w2 = waitbar(0, ' Creating points');
i1 = xpix.*0;

% create a vector of different Ns from 0 to nmax
ns = rand(round(m*0.25),1)*(nmax-nmin)/2 + nmin;
ss = rand(round(m*0.25),1)*(sigmax - sigmin)/2 + sigmin;
% Create a gaussian
% i1 = ntrue.*(2*pi*sigma^2)^-1.*exp(-((xpix-x0true).^2 +(ypix-y0true).^2)./(2*sigma.^2))+B;
%             i1x = 1/2.*(erf((xpix - x0true + 1/2)./(2*sigma^2)^0.5)-erf((xpix - x0true - 1/2)./(2*sigma^2)^0.5)); % error function of x integration over psf for each pixel
%             i1y = 1/2.*(erf((ypix - y0true + 1/2)./(2*sigma^2)^0.5)-erf((ypix - y0true - 1/2)./(2*sigma^2)^0.5)); % error function of y integration over psf for each pixel
%    i1 = i1x.*i1y;
for i = 1: round(m*0.25)
    sigma = ss(i);
    i1x = 1/2.*(erf((xpix - x0true + 1/2)./(2*sigma^2)^0.5)-erf((xpix - x0true - 1/2)./(2*sigma^2)^0.5)); % error function of x integration over psf for each pixel
    i1y = 1/2.*(erf((ypix - y0true + 1/2)./(2*sigma^2)^0.5)-erf((ypix - y0true - 1/2)./(2*sigma^2)^0.5)); % error function of y integration over psf for each pixel
    i1 = i1x.*i1y;
    i2(:,:) = double(imnoise(uint16(ns(i)*i1+B), 'poisson'))+.00001;
    x1 = [x1 ; i2(:).'];
    y1 = [y1; 1];
end