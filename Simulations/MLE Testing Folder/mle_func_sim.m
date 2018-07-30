clear all
close all
clc

%% Initialization
global xpix ypix wbox
rbox = 4;
[xpix, ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
% xpix = gpuArray(xpix);
% ypix = gpuArray(ypix);
wbox = 2*rbox+1;
i1 = xpix.*0;

% variables of image
B = 1;
ntrue = 100;
x0true = 2*(rand-0.5);
y0true = 2*(rand - 0.5);
% x0true = 0;
% y0true = 0;
sigma2 = 1.75;   %this would be rayleigh radius in pixel space
sigx = 1.18;   % sigma used for fitting
sigy = 1.24;
% r0t_pix = 2*sigma; % 1/e^2 radius
r0t_um = .2;        % Rayleigh radius in um
pix2pho = 1;
q = r0t_um / (sigma2);
wbox_um = q*wbox;

frames = 10000;
w2 = waitbar(0, ' Creating points');
i1 = xpix.*0;
% Create a gaussian
% i1 = ntrue.*(2*pi*sigma^2)^-1.*exp(-((xpix-x0true).^2 +(ypix-y0true).^2)./(2*sigma.^2))+B;
            i1x = 1/2.*(erf((xpix - x0true + 1/2)./(2*sigx^2)^0.5)-erf((xpix - x0true - 1/2)./(2*sigx^2)^0.5)); % error function of x integration over psf for each pixel
            i1y = 1/2.*(erf((ypix - y0true + 1/2)./(2*sigy^2)^0.5)-erf((ypix - y0true - 1/2)./(2*sigy^2)^0.5)); % error function of y integration over psf for each pixel
   i1 = ntrue * i1x.*i1y+B;                 

%% Create Frames with noise
for i = 1:frames
    waitbar(i/frames,w2, 'Creating points');
    i2(:,:,i) = single(imnoise(uint16(i1), 'poisson'))+.00001;
%     i2(:,:,i) = i1;
end
xf = [];
yf = [];
N = [];
O = [];
sx = [];
sy = [];

for i = 1:frames
    [fits] = func_mle_crlb(i2(:,:,i),0,0,3,0);
    xf = [xf;fits(1)];
    yf = [yf;fits(2)];
    N = [N;fits(3)];
    O = [O;fits(6)];
    sx = [sx;fits(4)];
    sy = [sy;fits(5)];
end

% remove 'failed' locs
ind = abs(xf) >= rbox | abs(yf) >= rbox | abs(sx) >= rbox | abs(sy) >= rbox |;
xf(ind) = [];
yf(ind) = [];
N(ind) = [];
O(ind) = [];
sx(ind) = [];
sy(ind) = [];

    