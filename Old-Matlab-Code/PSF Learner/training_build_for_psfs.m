%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comprehensive training set
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc


max_i = 10000;
max_j = 100;  % max number of noised images of a certain parameter set
[xpix, ypix]= meshgrid(-3:3,-3:3);
x0true = 0;
y0true = 0;
max_N0 = 300;
sigma2 = 2;   %this would be rayleigh radius in pixel space
sigma0 = sigma2/2;   % gaussian sigma
B0 = 3;
% w2 = waitbar(0, ' Creating points');
i1 = xpix.*0;
% Create a gaussian
y = -1*ones(max_i*max_j,1);
i3 = zeros(max_i*max_j,49);
% tic
for i = 1:max_i
    xtrue = 8*(rand-0.5);
    ytrue = 8*(rand-0.5);
    N = randi(max_N0);
    sigma = 0.05 * randn + sigma0;
    B = abs(randn + B0);
    i1 = N.*(2*pi*(sigma)^2)^-1.*exp(-((xpix-xtrue).^2 +(ypix-ytrue).^2)./(2*(sigma).^2))+B^2;
    for j = 1:max_j
        i2 = imnoise(uint8(i1),'poisson');
        i3((i-1)*max_j + j,:) = double(i2(:)).';
        if round(xtrue) == 0 && round(ytrue) == 0
            y((i-1)*max_j + j,1) = 1;
        else
            y((i-1)*max_j + j,1) = 0;
        end
%         tim = toc; %this is your tau
%         per = ((i-1)*max_j+j)/(max_i*max_j);
%         m = tim/per;
%         hour = floor(m/3600);
%         mins = floor((m-hour*3600)/60);
%         sec = round(m - hour*3600 - mins*60);
%         
%         disp([num2str(100*((i-1)*max_j+j)/(max_i*max_j)), '% complete estimated time ', num2str(mins),' minutes ', num2str(sec), ' seconds'])
    end
    disp(num2str(i/max_i));
end
save('Training_set.mat','i3','y');