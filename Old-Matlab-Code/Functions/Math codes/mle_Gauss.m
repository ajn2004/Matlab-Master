function [xf,yf,sx,sy,N,O] = mle_Gauss(i3)

[m,n,o] = size(i3);
if o > 1
    error('Only 1 Frame at a time');
elseif m ~=n
    error('square images only');
end
i3 = i3./max(i3(:));
[xpix, ypix] = meshgrid(-(n-1)/2:(n-1)/2,-(m-1)/2:(m-1)/2);
beta0 = [(sum(sum(xpix.*i3))/sum(i3(:))), (sum(sum(ypix.*i3))/sum(i3(:))), sum(i3(:)), 1.5, 1.5, min(i3(:))];
for k = 1:20

    Ex = 1/2.*(erf((xpix - beta0(1) + 1/2)./(2*beta0(4)^2)^0.5)-erf((xpix - beta0(1) - 1/2)./(2*beta0(4)^2)^0.5)); % error function of x integration over psf for each pixel
    Ey = 1/2.*(erf((ypix - beta0(2) + 1/2)./(2*beta0(5)^2)^0.5)-erf((ypix - beta0(2) - 1/2)./(2*beta0(5)^2)^0.5)); % error function of y integration over psf for each pixel
       
    u(:,:,k) = beta0(3).*Ex.*Ey + beta0(6);

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
    
    % update variables
    beta0(1) = beta0(1) - sum(sum(dudx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udx2.*((i3./u(:,:,k))-1) - dudx.^2.*i3./(u(:,:,k).^2))));
    beta0(2) = beta0(2) - sum(sum(dudy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udy2.*((i3./u(:,:,k))-1) - dudy.^2.*i3./(u(:,:,k).^2))));
    beta0(3) = beta0(3) - sum(sum(dudi.*((i3./u(:,:,k))-1)))/(sum(sum(d2udi2.*((i3./u(:,:,k))-1) - dudi.^2.*i3./(u(:,:,k).^2))));
    beta0(4) = beta0(4) - sum(sum(dudsx.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsx2.*((i3./u(:,:,k))-1) - dudsx.^2.*i3./(u(:,:,k).^2))));
    beta0(5) = beta0(5) - sum(sum(dudsy.*((i3./u(:,:,k))-1)))/(sum(sum(d2udsy2.*((i3./u(:,:,k))-1) - dudsy.^2.*i3./(u(:,:,k).^2))));
    beta0(6) = beta0(6) - sum(sum(dudb.*((i3./u(:,:,k))-1)))/(sum(sum(d2udb2.*((i3./u(:,:,k))-1) - dudb.^2.*i3./(u(:,:,k).^2))));

end
xf = beta0(1);
yf = beta0(2);
sx = beta0(3);
sy = beta0(4);
N  = beta0(5)*max(i3(:));
O  = beta0(6)*max(i3(:));