function e1 = func_gauss_blur(i1, ranged, type)


%Gaussian blur and edge detection

% clear all
% close all
% clc
%
% [fname, fpath] = uigetfile('*tif'); % choose an image
%
% i1 = rgb2gray(imread([fpath,fname])); % convert image to grayscale

fi1 = fftshift(fft2(i1));% get fft of image
%
[m, n] = size(i1);
[X, Y] = meshgrid(1:n,1:m);

for k = ranged % loop over fourier gaussian widths
    
    %Convolution of a gaussian is a guassian
    Z = exp(-2*((X-floor(n/2)-1).^2 + (Y -floor(m/2)-1).^2)./(k^2));
    ffi1 = fi1.*Z; %Convolutions become multiplications in fourier space
    
    %Radial Threshold on FFT
    %  k = 400;
    %  for i = 1:numel(fi1(:,1))
    %     for j = 1:numel(fi1(1,:))
    %         if ((i-floor(m/2)-1)^2 + (j-floor(n/2)-1)^2)^0.5 > k
    %             ffi1(i,j) = 0;
    %         end
    %     end
    %  end
    
    %Converted FFT
    i1f = abs(ifft2(fftshift(ffi1)));
    
    % Edge Filter
    e1 = edge(i1f,type);
    
    %Build Fame
    imagesc([i1, abs(ffi1); i1f, double(e1)*70],[0,70])
    colormap('gray')
    axis image
    title('Original Image                2D FFT .* Gaussian')
    xlabel('Reconstructed Image                Edge Filter Result');
    ylabel(['e^-2 radius = ', num2str(k),' pixels'])
    drawnow
    % Grab Frame
%     M(k) = getframe(gcf);
%     pause(0.2);
end
% Create Gif
% movie2gif(M,[type,'_edge.gif'],'DelayTime',0.2);