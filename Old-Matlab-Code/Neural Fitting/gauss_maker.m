%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Maker
%
% It Makes Gaussians and save the results
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% variables of image
global xpix ypix wbox
rbox = 4;
[xpix, ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
frames = 10000;
w2 = waitbar(0, ' Creating points');

for i = 1:frames
    waitbar(i/frames,w2, 'Creating points');
    B(i,1) = 3*rand;
    ntrue(i,1) = 500*rand+50;
    % x0true = 2*(rand-0.5);
    % y0true = 2*(rand - 0.5);
    x0true(i,1) = rand - 0.5;
    y0true(i,1) = rand - 0.5;
    sigma2(i,1) = 1.75+2*rand;   %this would be rayleigh radius in pixel space
%     sigr = 
    sigx(i,1) = (1)*sigma2(i)/2;   % sigma used for fitting
    sigy(i,1) = (1)*sigma2(i)/2;
    % Create a gaussian
    % i1 = ntrue.*(2*pi*sigma^2)^-1.*exp(-((xpix-x0true).^2 +(ypix-y0true).^2)./(2*sigma.^2))+B;
    i1x = 1/2.*(erf((xpix - x0true(i) + 1/2)./(2*sigx(i)^2)^0.5)-erf((xpix - x0true(i) - 1/2)./(2*sigx(i)^2)^0.5)); % error function of x integration over psf for each pixel
    i1y = 1/2.*(erf((ypix - y0true(i) + 1/2)./(2*sigy(i)^2)^0.5)-erf((ypix - y0true(i) - 1/2)./(2*sigy(i)^2)^0.5)); % error function of y integration over psf for each pixel
    i1 = ntrue(i) * i1x.*i1y+B(i);
    
    %% Create Frames with noise
    
    i3(:,:) = double(imnoise(uint16(i1), 'poisson'))+.00001;
    i2(i,:) = i3(:).';
    %     i2(:,:,i) = i1;
end
close(w2);
truth = [x0true, y0true, ntrue, sigx, sigy, B];
clearvars -except i2 truth
save('Test_gauss.mat');
