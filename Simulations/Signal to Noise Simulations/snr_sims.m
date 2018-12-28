close all; clearvars; clc; % clean up
% SNR_Localization_simulation
% This simulation will take fitted parameters of latex beads and use them
% to create simulated point spread functions which are used for subsequent
% fitting and comparison of localization information

%% User Settings
SNRs = 10.^(0:0.1:2); % these will enter the N to offset calculation via the equation snr = N/offset^0.5
frames = 1000; % Number of molecules to create at each setting
pixw = 4;   % size of total image molecule is created on
q = 0.128;
% End User Settings

%% PSF Developement
load('PSF_params.mat'); % Load bead data from a ladder test

% Choose fitting parameters


x0 = rand(1) - 0.5; % Center the Molecule somewhere in the first pixel
y0 = rand(1) - 0.5; % Center the Molecule somewhere in the first pixel
N = 2000; % SNR will be created from setting the offset relative to the N

% Allow user to 'choose' a sigma value based on height of the molecule
% fitted
figure('Units','Normalized','Outerposition', [0.25 0.25 0.5 0.5]) % Create a full frame figure
plot(coords(:,3)*q*1000); % Show stair step pattern 
xlabel('Framenumber'); % proper label of graph
ylabel('Axial position in pixels');% proper label of graph
title('Choose a height to Model'); % instruction to user
[x,y] = ginput(1); % Graphically chose the parameter desired
id = round(x); % Choose molecule index determined by user
sx = fits(id,4); % assign sigma values from molecule id
sy = fits(id,5);
B = 10;
[xpix,ypix] = meshgrid(-pixw:pixw); % Create mesh grids for desired window
for i = 1:numel(SNRs) % Loop over SNRs
    % Variables to Loop Over
%     i1x = 1/2.*(erf((xpix - x0 + 1/2)./(2*3.92^2)^0.5)-erf((xpix - x0 - 1/2)./(2*3.92^2)^0.5)); % error function of x integration over psf for each pixel
%     i1y = 1/2.*(erf((ypix - y0 + 1/2)./(2*3.92^2)^0.5)-erf((ypix - y0 - 1/2)./(2*3.92^2)^0.5)); % error function of y integration over psf for each pixel
id = randi(numel(fits(:,1)));
sx = fits(id,4); % assign sigma values from molecule id
sy = fits(id,5);
%     B = N/SNRs(i).*i1x.*i1y; % SNR = N/B algebra gives us the code
    B0{i} = B;
        i1 = xpix.*0;
        % Create a gaussian
        i1x = 1/2.*(erf((xpix - x0 + 1/2)./(2*sx^2)^0.5)-erf((xpix - x0 - 1/2)./(2*sx^2)^0.5)); % error function of x integration over psf for each pixel
        i1y = 1/2.*(erf((ypix - y0 + 1/2)./(2*sy^2)^0.5)-erf((ypix - y0 - 1/2)./(2*sy^2)^0.5)); % error function of y integration over psf for each pixel
        i1 = N * i1x.*i1y+B;
        im0(:,:,i) = i1;
    for j = 1:frames
        im1(:,:,j) = double(imnoise(uint16(i1),'poisson'));
%           im1(:,:,j) = double(imnoise(uint16(i1),'poisson'));
%           im1(:,:,j) = double(imnoise(uint16(i1),'poisson')) - imgaussfilt(B,1.5);
          imj = im1(:,:,j);
          sim = sort(imj(:));
%           im1(:,:,j) = im1(:,:,j) - mean(sim(1:4));
%         surf(im1(:,:,j))
%         drawnow
    end
    A{i} = im1.*(im1>0);
end


clear im1









load('z_calib.mat') % load latex calibrated data
cord0 = coords;


for i = 1:numel(A)
   iloc = A{i};
   [fits, crlbs, llv] = slim_locs(iloc);
   zf = getdz(fits(:,4),fits(:,5),cal.z_cal)/q;
   coords = [fits(:,1:2),zf];
   crl{i} = crlbs;
   cords{i} = coords;
   fit{i} = fits;
   try
       xstd(i) = std(coords(:,1)*q*1000);
       ystd(i) = std(coords(:,2)*q*1000);
       sxstd(i) = std(fits(:,4)*q);
       systd(i) = std(fits(:,5)*q);
       zstd(i) = std(coords(:,3)*q*1000);
   catch
   end
end     
plot(SNRs,zstd)
hold on
plot(SNRs,xstd)
plot(SNRs,ystd)
set(gca,'XScale','Log')
legend('Z-precision','X-Precision','Y-Precision')
xlabel('SNR')
ylabel('Localization Prescion in nm')
title('Subtracting Constant Background with noisy Guassian offset')

for i = 1:numel(crl)
    uncx(i) = mean(crl{i}(:,1).^0.5*128);
    uncy(i) = mean(crl{i}(:,2).^0.5*128);
    stx(i) = std(fit{i}(:,1)*128);
    sty(i) = std(fit{i}(:,2)*128);
end
% figure
% plot(SNRs,sxstd)
% hold on
% plot(SNRs,systd)
% legend('Sigx','Sigy')
% set(gca,'XScale','Log')