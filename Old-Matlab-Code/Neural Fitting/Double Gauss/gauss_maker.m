%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Maker
%
% It Makes Gaussians and save the results
% This script is adapted to create multiple gaussians
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% variables of image
global xpix ypix wbox
rbox = 4;
[xpix, ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
frames = 3000;
w2 = waitbar(0, ' Creating points');

for i = 1:frames
    if round(i/25) == i/25
        gnum(i,1) = 1;
    else
        gnum(i,1)=0;
    end
%     gnum(i,1) = randi(3) - 1;
    waitbar(i/frames,w2, 'Creating points');
    B(i,1) = 10*rand;
    i1 = zeros(rbox*2+1,rbox*2+1)+ B(i,1);
    ntrue(i,1) = 500*rand+50;
    % x0true = 2*(rand-0.5);
    % y0true = 2*(rand - 0.5);
    x0true(i,1) = 2*rand - 1;
    
    y0true(i,1) = 2*rand - 1;
    
    sigma0(i,1) = 1.75+2*rand;   %this would be rayleigh radius in pixel space
    %     sigr =
    sigx0(i,1) = (1)*sigma0(i)/2;   % sigma used for fitting
    sigy0(i,1) = (1)*sigma0(i)/2;
    x1true(i,1) = 2*rand - 1;
    y1true(i,1) = 2*rand - 1;
    sigma1(i,1) = 1.75+2*rand;   %this would be rayleigh radius in pixel space
    %     sigr =
    sigx1(i,1) = (1)*sigma1(i)/2;   % sigma used for fitting
    sigy1(i,1) = (1)*sigma1(i)/2;
    if gnum(i) > 0
        
        % Create a gaussian
        % i1 = ntrue.*(2*pi*sigma^2)^-1.*exp(-((xpix-x0true).^2 +(ypix-y0true).^2)./(2*sigma.^2))+B;
        i1x = 1/2.*(erf((xpix - x0true(i) + 1/2)./(2*sigx0(i)^2)^0.5)-erf((xpix - x0true(i) - 1/2)./(2*sigx0(i)^2)^0.5)); % error function of x integration over psf for each pixel
        i1y = 1/2.*(erf((ypix - y0true(i) + 1/2)./(2*sigy0(i)^2)^0.5)-erf((ypix - y0true(i) - 1/2)./(2*sigy0(i)^2)^0.5)); % error function of y integration over psf for each pixel
        i1 = i1 + ntrue(i) * i1x.*i1y;
        
        % Create second gaussian
        if gnum(i) == 2
            
            i2x = 1/2.*(erf((xpix - x1true(i) + 1/2)./(2*sigx1(i)^2)^0.5)-erf((xpix - x1true(i) - 1/2)./(2*sigx1(i)^2)^0.5)); % error function of x integration over psf for each pixel
            i2y = 1/2.*(erf((ypix - y1true(i) + 1/2)./(2*sigy1(i)^2)^0.5)-erf((ypix - y1true(i) - 1/2)./(2*sigy1(i)^2)^0.5)); % error function of y integration over psf for each pixel
            i1 = i1 + ntrue(i) * i2x.*i2y;
        end
    end
    %% Create Frames with noise
    
    i3(:,:) = double(imnoise(uint16(i1), 'poisson'))+.00001;
    i2(i,:) = i3(:).';
    %     i2(:,:,i) = i1;
    %     imagesc(i1);
    %     colormap('gray');
    %     hold on
    %     if gnum(i) > 0
    %         plot(x0true(i,1)+rbox+1,y0true(i,1)+rbox+1,'.r');
    %
    %         if gnum(i) == 2
    %             plot(x1true(i,1)+rbox+1,y1true(i,1)+rbox+1,'.b');
    %         end
    %
    %     end
    %     hold off
    %     drawnow
end
close(w2);
% truth = gnum;
truth = [gnum, x0true, y0true, x1true, y1true, ntrue, sigx0, sigy0, sigx1,sigy1 B];
clearvars -except i2 truth
save('Test_gauss.mat');
