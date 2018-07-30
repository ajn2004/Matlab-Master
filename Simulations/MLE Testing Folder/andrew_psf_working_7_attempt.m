%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andrew's Gaussian MLE Localization Code
% AJN 3/30/15
%
% This script will generate a stack of gaussian point spread functions with
% the variables B, ntrue, x0true, y0true, r0t_pix
% The outputs of this script are xf, yf, a0, off0 for x position, y
% position number of photons and offset (background) respectively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Initialization
global xpix ypix wbox
rbox = 4;
[xpix, ypix] = meshgrid(-rbox:rbox,-rbox:rbox);

wbox = 2*rbox+1;
i1 = xpix.*0;

% variables of image
B = 1;
ntrue = 100;
x0true = 2*(rand-0.5);
y0true = 2*(rand - 0.5);

sigma2 = 1.75;   %this would be rayleigh radius in pixel space
sigx = 1.18;   % sigma used for fitting
sigy = 1.24;

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
% imagesc(i2);
% colormap(gray);
% gpui2 = gpuArray(i2);

close (w2)
w1 = waitbar(0,'Fitting, be patient');

%% Find molecules and fit them
for l = 1:frames
    waitbar(l/frames,w1 , 'Fitting be patient')
% i3 = gpui2(:,:,l);
i3 = i2(:,:,l);
% i3 = i1;
[highrow, highcol] = find(i3 == max(i3(:)),1);
highpix = max(i3(:));
k=0;
zlin = i3(:);
xfake = zlin.*0;
clear xguess yguess peakguess r0_allguess offguess beta0 Ex Ey u dudx dudy dudi dudb d2udx2 d2udy2 d2udi2 d2udb2 lowx lowy hix hiy
xguess = xpix(highrow,highcol);
yguess = ypix(highrow,highcol);


lowx = round(xguess-2*sigma2+rbox+1);
hix = round(xguess+2*sigma2+rbox+1);
lowy = round(yguess-2*sigma2+rbox+1);
hiy = round(xguess+2*sigma2+rbox+1);


if round(xguess-2*sigma2+rbox+1) <=0
    lowx = 1;
end
if round(xguess+2*sigma2+rbox+1) >= max(max(xpix))
    hix = max(max(xpix));
end
if round(yguess-2*sigma2+rbox+1) <=0
    lowy = 1;
    end
if round(xguess+2*sigma2+rbox+1) >= max(max(ypix))
    hiy = max(max(ypix));
end
peakguess = max(max(i3))/2*pi*(2*sigma2/2)^2;
offguess = 0;

beta0 = [ xguess, yguess, peakguess, sigma2/2, sigma2/2, offguess];
% 
% [beta,R,J,CovB,MSE] = nlinfit(xfake,zlin,@gaussguess, beta0);
% elapsed_nlon = toc;
%% MLE approximation
sigma = beta0(4);

fittime(l) = 1;
% u = gpuArray(zeros(wbox,wbox,frames));
% while fittime(l) < 100

for k = 1:10
%     if beta0(5) <= 0
%         beta0(5) =0;
%     end
%Define psf and error function for each pixel
%     Ex = 1/2.*(erf((xpix - beta0(1) + 1/2)./(2*sigma^2)^0.5)-erf((xpix - beta0(1) - 1/2)./(2*sigma^2)^0.5)); % error function of x integration over psf for each pixel
%     Ey = 1/2.*(erf((ypix - beta0(2) + 1/2)./(2*sigma^2)^0.5)-erf((ypix - beta0(2) - 1/2)./(2*sigma^2)^0.5)); % error function of y integration over psf for each pixel
    
    Ex = 1/2.*(erf((xpix - beta0(1) + 1/2)./(2*beta0(4)^2)^0.5)-erf((xpix - beta0(1) - 1/2)./(2*beta0(4)^2)^0.5)); % error function of x integration over psf for each pixel
    Ey = 1/2.*(erf((ypix - beta0(2) + 1/2)./(2*beta0(5)^2)^0.5)-erf((ypix - beta0(2) - 1/2)./(2*beta0(5)^2)^0.5)); % error function of y integration over psf for each pixel

%    Ex = 1/2*(erf((xpix - beta0(1) + 1/2)/(2*sigma^2))-erf((xpix - beta0(1) - 1/2)/(2*sigma^2))); % error function of x integration over psf for each pixel
%    Ey = 1/2*(erf((ypix - beta0(2) + 1/2)/(2*sigma^2))-erf((ypix - beta0(2) - 1/2)/(2*sigma^2))); % error function of y integration over psf for each pixel


    u(:,:,k) = beta0(3).*Ex.*Ey + beta0(6);
%     surf(U,:,:
    % partial derivatives of variables of interest
    dudx = beta0(3)*(2*pi*beta0(4)^2)^-0.5.*(exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1)-exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey;
    dudy = beta0(3)*(2*pi*beta0(5)^2)^-0.5.*(exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1)-exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex;
    dudsx = beta0(3)*(2*pi)^(-1/2)*beta0(4)^(-2).*((xpix - beta0(1) - 1/2).*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey; % pd sigx
    dudsy = beta0(3)*(2*pi)^(-1/2)*beta0(5)^(-2).*((ypix - beta0(2) - 1/2).*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2) .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex; % pd sigy
    dudi = Ex.*Ey;
    dudb = 1;
    
    % Second partial derivatives of variables of interest
    d2udx2 = beta0(3)*(2*pi)^-0.5*beta0(4)^-3*((xpix - beta0(1) - 1/2).*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey;
    d2udy2 = beta0(3)*(2*pi)^-0.5*beta0(5)^-3*((ypix - beta0(2) - 1/2).*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2) .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex;
        d2udsx2 = beta0(3).*Ey.*(2*pi)^-0.5.*((beta0(4)^-5.* ((xpix - beta0(1) - 1/2).^3.*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2).^3.*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1))) ...
            - 2.*beta0(4).^-3.*((xpix - beta0(1) - 1/2).*   exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*  exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)));
        % second partial for sigmay
        d2udsy2 = beta0(3).*Ex.*(2*pi)^-0.5.*((beta0(5)^-5.* ((ypix - beta0(2) - 1/2).^3.*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2).^3.*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1))) ...
            - 2.*beta0(5).^-3.*((ypix - beta0(2) - 1/2).*   exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2)   .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)));

    
    d2udi2 = 0;
    d2udb2 = 0;
    
    xlap(k,l) = beta0(1);
    ylap(k,l) = beta0(2);
    siglap(k,l) = beta0(4);
    nlap(k,l) = beta0(3);
    bglap(k,l) = beta0(6);
    % update variables
    beta0(1) = beta0(1) - sum(sum(dudx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udx2.*((i3./u(:,:,k))-1) - dudx.^2.*i3./(u(:,:,k).^2))));
    beta0(2) = beta0(2) - sum(sum(dudy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udy2.*((i3./u(:,:,k))-1) - dudy.^2.*i3./(u(:,:,k).^2))));
    beta0(3) = beta0(3) - sum(sum(dudi.*((i3./u(:,:,k))-1)))/(sum(sum(d2udi2.*((i3./u(:,:,k))-1) - dudi.^2.*i3./(u(:,:,k).^2))));
    beta0(4) = beta0(4) - sum(sum(dudsx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsx2.*((i3./u(:,:,k))-1) - dudsx.^2.*i3./(u(:,:,k).^2))));
    beta0(5) = beta0(5) - sum(sum(dudsy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsy2.*((i3./u(:,:,k))-1) - dudsy.^2.*i3./(u(:,:,k).^2))));
    beta0(6) = beta0(6) - sum(sum(dudb.*((i3./u(:,:,k))-1)))/(sum(sum(d2udb2.*((i3./u(:,:,k))-1) - dudb.^2.*i3./(u(:,:,k).^2))));
    


end   
    xlap(k+1,l) = beta0(1);
    ylap(k+1,l) = beta0(2);
    siglap(k+1,l) = beta0(4);
    nlap(k+1,l) = beta0(3);
    bglap(k+1,l) = beta0(6);
    p(l) = sum(sum(i3.*real(log(u(:,:,end)))-u(:,:,end)-i3.*real(log(i3))+i3));

% end

xf(l) = beta0(1);
yf(l) = beta0(2);
a0(l) = beta0(3);
r0(l) = (beta0(4)*beta0(5))^0.5;
sigx(l) = beta0(4);
sigy(l) = beta0(5);
off0(l) = beta0(6);
end

close (w1)
% xf = beta0(1);
% yf = beta0(2);
% a0 = beta0(3);
% r0 = beta0(4);
% off0 = beta0(5);
% rf = (xf.^2+yf.^2).^0.5;
% hist(rf*q*1000,1:2:80);
% drawnow
% % end
% imagesc(u*1000);
% subplot(1,2,1);imagesc(i2);
% subplot(1,2,2);imagesc(u(:,:,l)*1000000);
% axis image
% colormap('gray');
% plot(1:60,xf.*q)
numel(find(fittime==101))
imagesc(i2(:,:,1));
axis image
colormap('gray')
hold on
plot(xf+rbox+1, yf+rbox+1, '.b');
plot(x0true+rbox+1,y0true+rbox+1,'.r')
title('PSF with blue localizations and red source')
figure
plot(1:k+1,siglap-1.18);
title('Convergence of fitting parameter')
xlabel('Iteration')
ylabel('Value')
figure
histogram(xf-x0true);