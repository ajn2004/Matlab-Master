%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Molecule frame simulator
% Simulates a frame of molecules
%
%
% AJN 10/23/15
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

num_mol = randi(3);
max_pho = 200;
bkg = 1;

pixx = 10; % number of pixels in x
pixy = 10; %number of pixels in y
q = 133; % pixel size in nanometers
NA = 1.45;  %numerical aperture of lens
lam = 573; %wavelength of peak emission in nanometers

sigma2 = 0.55*lam/NA; % 1/e^2 radius of molecular image
sigma = sigma2/2;
sig2pix = sigma2/q;  % radius in pixels
sigpix = sig2pix/2;

i1 = zeros(q*pixy,q*pixx);

for i = 1:num_mol
    index(i,1) = randi(q*pixx);
    index(i,2) = randi(q*pixy);
    N(i) = randi(max_pho);
    i1(index(i,2),index(i,1)) = N(i);
end

% imagesc(i1);
% axis image

%% build gaussian for convolution
[X, Y] = meshgrid(-3*q:3*q,-3*q:3*q);
gauss = (2*pi*sigma)^-2.*exp(-(X.^2 + Y.^2)./(2*sigma^2));
i2 = conv2(i1,gauss,'same');


i3 = zeros(pixy,pixx);
pos = zeros(pixy, pixx);
%% build pixel frame
for i = 1:pixx
    for j = 1:pixy
       i3(j,i) = sum(sum(i2((j-1)*q+1:j*q,(i-1)*q+1:i*q)));
       
    end
end

%% Add noise in
i3 = i3 + bkg;
i4 = imnoise(uint16(i3),'poisson');

figure
title('x marks pixel that containes the molecule');
subplot(2,2,1);
imagesc(i4)
title('x marks pixel that contains the molecule');
axis image
colormap(gray)
hold on
plot(ceil(index(:,1)/q),ceil(index(:,2)/q),'.r','MarkerSize',20,'Marker','x')
subplot(2,2,2);
imagesc(i3);
title(['Image of ', num2str(num_mol),' molecules']);
axis image
colormap(gray)
hold on
plot(ceil(index(:,1)/q),ceil(index(:,2)/q),'.r','MarkerSize',20,'Marker','x')
subplot(2,2,3:4)
imagesc(i4);
axis image
colormap(gray)