function [xf, yf,N, sigx, sigy, offset] = func_mle(i3, xguess, yguess, sigma2)
% Function to perform the MLE calculation on an image
[xpix, ypix] = meshgrid(-(numel(i3(:,1))-1)/2:(numel(i3(:,1))-1)/2,-(numel(i3(:,1))-1)/2:(numel(i3(:,1))-1)/2 );
peakguess = max(max(i3))/2*pi*(2*1.5)^2;
offguess = min(i3(:));

beta0 = [ xguess, yguess, peakguess, sigma2, sigma2, offguess];
% 
% [beta,R,J,CovB,MSE] = nlinfit(xfake,zlin,@gaussguess, beta0);
% elapsed_nlon = toc;
%% MLE approximation
sigma = beta0(4);


% u = gpuArray(zeros(wbox,wbox,frames));
% while fittime(l) < 100

for k = 1:45
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
    dudx = beta0(3)*(2*pi*sigma^2)^-0.5.*(exp(-(xpix -beta0(1) - 1/2).^2.*(2*sigma^2)^-1)-exp(-(xpix -beta0(1) + 1/2).^2.*(2*sigma^2)^-1)).*Ey;
    dudy = beta0(3)*(2*pi*sigma^2)^-0.5.*(exp(-(ypix -beta0(2) - 1/2).^2.*(2*sigma^2)^-1)-exp(-(ypix -beta0(2) + 1/2).^2.*(2*sigma^2)^-1)).*Ex;
    dudsx = beta0(3)*(2*pi)^(-1/2)*beta0(4)^(-2).*((xpix - beta0(1) - 1/2).*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey; % pd sigx
    dudsy = beta0(3)*(2*pi)^(-1/2)*beta0(5)^(-2).*((ypix - beta0(2) - 1/2).*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2) .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex; % pd sigy
    dudi = Ex.*Ey;
    dudb = 1;
    
    % Second partial derivatives of variables of interest
    d2udx2 = beta0(3)*(2*pi)^-0.5*sigma^-3*((xpix - beta0(1) - 1/2).*exp(-(xpix -beta0(1) - 1/2).^2.*(2*sigma^2)^-1) - (xpix - beta0(1) + 1/2) .*exp(-(xpix -beta0(1) + 1/2).^2.*(2*sigma^2)^-1)).*Ey;
    d2udy2 = beta0(3)*(2*pi)^-0.5*sigma^-3*((ypix - beta0(2) - 1/2).*exp(-(ypix -beta0(2) - 1/2).^2.*(2*sigma^2)^-1) - (ypix - beta0(2) + 1/2) .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*sigma^2)^-1)).*Ex;
        d2udsx2 = beta0(3).*Ey.*(2*pi)^-0.5.*((beta0(4)^-5.* ((xpix - beta0(1) - 1/2).^3.*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2).^3.*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1))) ...
            - 2.*beta0(4).^-3.*((xpix - beta0(1) - 1/2).*   exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*  exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)));
        % second partial for sigmay
        d2udsy2 = beta0(3).*Ex.*(2*pi)^-0.5.*((beta0(5)^-5.* ((ypix - beta0(2) - 1/2).^3.*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2).^3.*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1))) ...
            - 2.*beta0(5).^-3.*((ypix - beta0(2) - 1/2).*   exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2)   .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)));

    
    d2udi2 = 0;
    d2udb2 = 0;
    
    % update variables
    beta0(1) = beta0(1) - sum(sum(dudx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udx2.*((i3./u(:,:,k))-1) - dudx.^2.*i3./(u(:,:,k).^2))));
    beta0(2) = beta0(2) - sum(sum(dudy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udy2.*((i3./u(:,:,k))-1) - dudy.^2.*i3./(u(:,:,k).^2))));
    beta0(3) = beta0(3) - sum(sum(dudi.*((i3./u(:,:,k))-1)))/(sum(sum(d2udi2.*((i3./u(:,:,k))-1) - dudi.^2.*i3./(u(:,:,k).^2))));
    beta0(4) = beta0(4) - sum(sum(dudsx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsx2.*((i3./u(:,:,k))-1) - dudsx.^2.*i3./(u(:,:,k).^2))));
    beta0(5) = beta0(5) - sum(sum(dudsy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsy2.*((i3./u(:,:,k))-1) - dudsy.^2.*i3./(u(:,:,k).^2))));
    beta0(6) = beta0(6) - sum(sum(dudb.*((i3./u(:,:,k))-1)))/(sum(sum(d2udb2.*((i3./u(:,:,k))-1) - dudb.^2.*i3./(u(:,:,k).^2))));
    


end   
%     xlap(k+1,l) = beta0(1);
%     ylap(k+1,l) = beta0(2);
%     siglap(k+1,l) = beta0(4);
%     nlap(k+1,l) = beta0(3);
%     bglap(k+1,l) = beta0(6);
%     p(l) = sum(sum(i3.*real(log(u(:,:,end)))-u(:,:,end)-i3.*real(log(i3))+i3));

% end

xf = beta0(1);
yf = beta0(2);
N = beta0(3);
sigx = beta0(4);
sigy = beta0(5);
offset =  beta0(6);
end