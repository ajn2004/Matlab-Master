%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circular Density Tester
%
% This script will generate several circles randomly filled in and not
% filled in. The user will be able to provide several parameters such as
% size of circle, thickness of edge, fill rate, and background. The process
% will be done by drawing a circle with 1 pixel width, then convolving with
% a guassian. The main goal of this script is to develop the mathematical
% framework necessary to perform the radial density measurement
%
% THIS DOES NOT CONTAIN GUSTO!
% 
% AJN 2/19/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
figure('units','normalized','outerposition',[0 0 1 1])
%% Parameters to control
fact = 3;
% circle creation
for radius = 0:40
% radius = randi(20); % radius in pixels
xc = 0; % center of the circle
yc = 0; % center of circle

% Gaussian control
r0 = 5;  % radius of gaussian
a0 = 100; % peak value of gaussian
n = 2*pi*(r0/2)^2;

% Image parameters
back = 4; % background
pixs = 2000; % length of image to be created (if this is even the image length will be 1+ this value)

%% Image creation

% determine size of image
if pixs / 2 == round(pixs/2)
    sizes = -pixs/2 : pixs/2;
else
    pixe = pixs - 1;
    sizes = -pixe/2:pixe/2;
end

[X, Y] = meshgrid(sizes,sizes);
i1 = zeros(size(X));

im1 = circle_maker(xc,yc, radius,i1); % create logical image of a circle

% imagesc(im1); % show circle

gauss = a0/n * exp(-2*((X.*X)+(Y.*Y))/r0^2); %create a guassian whose area is a0

% figure

% convolve circle with gaussian
im2 = round(conv2(im1,gauss,'same'));

BW2 = logical(round(im2));
s = regionprops(BW2,'centroid','EquivDiameter');
% imagesc(im2);

%% This section is what measures the radial density function
cents = s.Centroid;
radius1 = s.EquivDiameter/2;

[rs, ds] = circular_density(cents(1),cents(2),radius1,im2,fact);
% figure
subplot(1,2,1);plot(rs/radius1,ds/max(ds(ds<10000000))); title('radial density');xlabel('radius in pixels');ylabel('Density');
xlim([0 3])
subplot(1,2,2); imagesc(im2);
drawnow
% pause(0.1)
end