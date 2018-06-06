function [fits, crlbs, llv] = func_mle_crlb(i3, xguess, yguess, sigma2)
% Function to perform the MLE calculation on an image
[xpix, ypix] = meshgrid(-(numel(i3(:,1))-1)/2:(numel(i3(:,1))-1)/2,-(numel(i3(:,1))-1)/2:(numel(i3(:,1))-1)/2 );
peakguess = max(max(i3))/2*pi*(2*1.5)^2;
offguess = min(i3(:));

beta0 = [ xguess, yguess, peakguess, sigma2, sigma2, offguess];
%
% (beta,R,J,CovB,MSE) = nlinfit(xfake,zlin,@gaussguess, beta0);
% elapsed_nlon = toc;
%% MLE approximation
sigma = beta0(4);


% u = gpuArray(zeros(wbox,wbox,frames));
% while fittime(l) < 100

for k = 1:45
    
    %Define psf and error function for each pixel
    
    Ex = 1/2.*(erf((xpix - beta0(1) + 1/2)./(2*beta0(4)^2)^0.5)-erf((xpix - beta0(1) - 1/2)./(2*beta0(4)^2)^0.5)); % error function of x integration over psf for each pixel
    Ey = 1/2.*(erf((ypix - beta0(2) + 1/2)./(2*beta0(5)^2)^0.5)-erf((ypix - beta0(2) - 1/2)./(2*beta0(5)^2)^0.5)); % error function of y integration over psf for each pixel
    
    
    u(:,:) = beta0(3).*Ex.*Ey + beta0(6);
    
    % partial derivatives of variables of interest
    dudx = beta0(3)*(2*pi*beta0(4)^2)^-0.5.*(exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1)-exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey;
    dudy = beta0(3)*(2*pi*beta0(5)^2)^-0.5.*(exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1)-exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex;
    dudsx = beta0(3)*(2*pi)^(-1/2)*beta0(4)^(-2).*((xpix - beta0(1) - 1/2).*exp(-(xpix -beta0(1) - 1/2).^2.*(2*beta0(4)^2)^-1) - (xpix - beta0(1) + 1/2) .*exp(-(xpix -beta0(1) + 1/2).^2.*(2*beta0(4)^2)^-1)).*Ey; % pd sigx
    dudsy = beta0(3)*(2*pi)^(-1/2)*beta0(5)^(-2).*((ypix - beta0(2) - 1/2).*exp(-(ypix -beta0(2) - 1/2).^2.*(2*beta0(5)^2)^-1) - (ypix - beta0(2) + 1/2) .*exp(-(ypix -beta0(2) + 1/2).^2.*(2*beta0(5)^2)^-1)).*Ex; % pd sigy
    dudn = Ex.*Ey;
    dudo = 1;
    
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
    beta0(1) = beta0(1) - sum(sum(dudx.*((i3./u(:,:))-1)))/(sum(sum(d2udx2.*((i3./u(:,:))-1) - dudx.^2.*i3./(u(:,:).^2))));
    beta0(2) = beta0(2) - sum(sum(dudy.*((i3./u(:,:))-1)))/(sum(sum(d2udy2.*((i3./u(:,:))-1) - dudy.^2.*i3./(u(:,:).^2))));
    beta0(3) = beta0(3) - sum(sum(dudn.*((i3./u(:,:))-1)))/(sum(sum(d2udi2.*((i3./u(:,:))-1) - dudn.^2.*i3./(u(:,:).^2))));
    beta0(4) = beta0(4) - sum(sum(dudsx.*((i3./u(:,:))-1)))/(sum(sum(d2udsx2.*((i3./u(:,:))-1) - dudsx.^2.*i3./(u(:,:).^2))));
    beta0(5) = beta0(5) - sum(sum(dudsy.*((i3./u(:,:))-1)))/(sum(sum(d2udsy2.*((i3./u(:,:))-1) - dudsy.^2.*i3./(u(:,:).^2))));
    beta0(6) = beta0(6) - sum(sum(dudo.*((i3./u(:,:))-1)))/(sum(sum(d2udb2.*((i3./u(:,:))-1) - dudo.^2.*i3./(u(:,:).^2))));
    
    if k ==45
        fisher(1) = sum(sum(dudx.*dudx ./ u));
        fisher(2) = sum(sum(dudx.*dudy ./ u));
        fisher(3) = sum(sum(dudx.*dudn ./ u));
        fisher(4) = sum(sum(dudx.*dudo ./ u));
        fisher(5) = sum(sum(dudx.*dudsx ./ u));
        fisher(6) = sum(sum(dudx.*dudsy ./ u));
        
        fisher(7) = sum(sum(dudy.*dudx ./ u));
        fisher(8) = sum(sum(dudy.*dudy ./ u));
        fisher(9) = sum(sum(dudy.*dudn ./ u));
        fisher(10) = sum(sum(dudy.*dudo ./ u));
        fisher(11) = sum(sum(dudy.*dudsx ./ u));
        fisher(12) = sum(sum(dudy.*dudsy ./ u));
        
        fisher(13) = sum(sum(dudn.*dudx ./ u));
        fisher(14) = sum(sum(dudn.*dudy ./ u));
        fisher(15) = sum(sum(dudn.*dudn ./ u));
        fisher(16) = sum(sum(dudn.*dudo ./ u));
        fisher(17) = sum(sum(dudn.*dudsx ./ u));
        fisher(18) = sum(sum(dudn.*dudsy ./ u));
        
        fisher(19) = sum(sum(dudo.*dudx ./ u));
        fisher(20) = sum(sum(dudo.*dudy ./ u));
        fisher(21) = sum(sum(dudo.*dudn ./ u));
        fisher(22) = sum(sum(dudo.*dudo ./ u));
        fisher(23) = sum(sum(dudo.*dudsx ./ u));
        fisher(24) = sum(sum(dudo.*dudsy ./ u));
        
        fisher(25) = sum(sum(dudsx.*dudx ./ u));
        fisher(26) = sum(sum(dudsx.*dudy ./ u));
        fisher(27) = sum(sum(dudsx.*dudn ./ u));
        fisher(28) = sum(sum(dudsx.*dudo ./ u));
        fisher(29) = sum(sum(dudsx.*dudsx ./ u));
        fisher(30) = sum(sum(dudsx.*dudsy ./ u));
        
        fisher(31) = sum(sum(dudsy.*dudx ./ u));
        fisher(32) = sum(sum(dudsy.*dudy ./ u));
        fisher(33) = sum(sum(dudsy.*dudn ./ u));
        fisher(34) = sum(sum(dudsy.*dudo ./ u));
        fisher(35) = sum(sum(dudsy.*dudsx ./ u));
        fisher(36) = sum(sum(dudsy.*dudsy ./ u));
        llv =  sum(sum(i3.*log(u) - u - i3.*log(i3) + i3));
    end
end
fisher = reshape(fisher,6,6);
det_fish = det(fisher);
xf_crlb = (fisher(8) * (fisher(15) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(21) * (fisher(16) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) + fisher(27) * (fisher(16) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(33) * (fisher(16) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) + fisher(28) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)))) - fisher(14) * (fisher(9) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(21) * (fisher(10) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) + fisher(27) * (fisher(10) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) - fisher(33) * (fisher(10) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)))) + fisher(20) * (fisher(9) * (fisher(16) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) - fisher(15) * (fisher(10) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) + fisher(27) * (fisher(10) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(16) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(33) * (fisher(10) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(16) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(18) - fisher(17) * fisher(12)))) - fisher(26) * (fisher(9) * (fisher(16) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(15) * (fisher(10) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) + fisher(21) * (fisher(10) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(16) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(33) * (fisher(10) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)) - fisher(16) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)) + fisher(22) * (fisher(11) * fisher(18) - fisher(17) * fisher(12)))) + fisher(32) * (fisher(9) * (fisher(16) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) + fisher(28) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(15) * (fisher(10) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) + fisher(21) * (fisher(10) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(16) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(27) * (fisher(10) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)) - fisher(16) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)) + fisher(22) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))))) / det_fish;
yf_crlb = -(-(fisher(1) * (fisher(15) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(21) * (fisher(16) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) + fisher(27) * (fisher(16) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(33) * (fisher(16) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) + fisher(28) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)))) - fisher(13) * (fisher(3) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(21) * (fisher(4) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(27) * (fisher(4) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) - fisher(33) * (fisher(4) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)))) + fisher(19) * (fisher(3) * (fisher(16) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) - fisher(15) * (fisher(4) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(27) * (fisher(4) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(16) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) - fisher(33) * (fisher(4) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(16) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(18) - fisher(17) * fisher(6)))) - fisher(25) * (fisher(3) * (fisher(16) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(34) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(15) * (fisher(4) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) + fisher(21) * (fisher(4) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(16) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) - fisher(33) * (fisher(4) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)) - fisher(16) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)) + fisher(22) * (fisher(5) * fisher(18) - fisher(17) * fisher(6)))) + fisher(31) * (fisher(3) * (fisher(16) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) + fisher(28) * (fisher(17) * fisher(24) - fisher(23) * fisher(18))) - fisher(15) * (fisher(4) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) + fisher(21) * (fisher(4) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(16) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) - fisher(27) * (fisher(4) * (fisher(17) * fisher(24) - fisher(23) * fisher(18)) - fisher(16) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)) + fisher(22) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))))) / det_fish);
N_crlb = (fisher(1) * (fisher(8) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(20) * (fisher(10) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) + fisher(26) * (fisher(10) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) - fisher(32) * (fisher(10) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)))) - fisher(7) * (fisher(2) * (fisher(22) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) + fisher(34) * (fisher(23) * fisher(30) - fisher(29) * fisher(24))) - fisher(20) * (fisher(4) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(26) * (fisher(4) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) - fisher(32) * (fisher(4) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)))) + fisher(19) * (fisher(2) * (fisher(10) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) - fisher(8) * (fisher(4) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(28) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(26) * (fisher(4) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) - fisher(10) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(32) * (fisher(4) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) - fisher(10) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(12) - fisher(11) * fisher(6)))) - fisher(25) * (fisher(2) * (fisher(10) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(34) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) - fisher(8) * (fisher(4) * (fisher(23) * fisher(36) - fisher(35) * fisher(24)) - fisher(22) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) + fisher(20) * (fisher(4) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) - fisher(10) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(34) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(32) * (fisher(4) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)) - fisher(10) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)) + fisher(22) * (fisher(5) * fisher(12) - fisher(11) * fisher(6)))) + fisher(31) * (fisher(2) * (fisher(10) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(28) * (fisher(11) * fisher(24) - fisher(23) * fisher(12))) - fisher(8) * (fisher(4) * (fisher(23) * fisher(30) - fisher(29) * fisher(24)) - fisher(22) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(24) - fisher(23) * fisher(6))) + fisher(20) * (fisher(4) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) - fisher(10) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(28) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(26) * (fisher(4) * (fisher(11) * fisher(24) - fisher(23) * fisher(12)) - fisher(10) * (fisher(5) * fisher(24) - fisher(23) * fisher(6)) + fisher(22) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))))) / det_fish;
off_crlb = -(-(fisher(1) * (fisher(8) * (fisher(15) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(33) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) - fisher(14) * (fisher(9) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(33) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) + fisher(26) * (fisher(9) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(15) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(33) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(32) * (fisher(9) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(15) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(27) * (fisher(11) * fisher(18) - fisher(17) * fisher(12)))) - fisher(7) * (fisher(2) * (fisher(15) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) + fisher(33) * (fisher(17) * fisher(30) - fisher(29) * fisher(18))) - fisher(14) * (fisher(3) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(26) * (fisher(3) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(15) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) - fisher(32) * (fisher(3) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(15) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(27) * (fisher(5) * fisher(18) - fisher(17) * fisher(6)))) + fisher(13) * (fisher(2) * (fisher(9) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(33) * (fisher(11) * fisher(30) - fisher(29) * fisher(12))) - fisher(8) * (fisher(3) * (fisher(29) * fisher(36) - fisher(35) * fisher(30)) - fisher(27) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(30) - fisher(29) * fisher(6))) + fisher(26) * (fisher(3) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) - fisher(9) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(32) * (fisher(3) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) - fisher(9) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(27) * (fisher(5) * fisher(12) - fisher(11) * fisher(6)))) - fisher(25) * (fisher(2) * (fisher(9) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(15) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) + fisher(33) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(8) * (fisher(3) * (fisher(17) * fisher(36) - fisher(35) * fisher(18)) - fisher(15) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) + fisher(14) * (fisher(3) * (fisher(11) * fisher(36) - fisher(35) * fisher(12)) - fisher(9) * (fisher(5) * fisher(36) - fisher(35) * fisher(6)) + fisher(33) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(32) * (fisher(3) * (fisher(11) * fisher(18) - fisher(17) * fisher(12)) - fisher(9) * (fisher(5) * fisher(18) - fisher(17) * fisher(6)) + fisher(15) * (fisher(5) * fisher(12) - fisher(11) * fisher(6)))) + fisher(31) * (fisher(2) * (fisher(9) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(15) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) + fisher(27) * (fisher(11) * fisher(18) - fisher(17) * fisher(12))) - fisher(8) * (fisher(3) * (fisher(17) * fisher(30) - fisher(29) * fisher(18)) - fisher(15) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(27) * (fisher(5) * fisher(18) - fisher(17) * fisher(6))) + fisher(14) * (fisher(3) * (fisher(11) * fisher(30) - fisher(29) * fisher(12)) - fisher(9) * (fisher(5) * fisher(30) - fisher(29) * fisher(6)) + fisher(27) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))) - fisher(26) * (fisher(3) * (fisher(11) * fisher(18) - fisher(17) * fisher(12)) - fisher(9) * (fisher(5) * fisher(18) - fisher(17) * fisher(6)) + fisher(15) * (fisher(5) * fisher(12) - fisher(11) * fisher(6))))) / det_fish);
sigx_crlb = (fisher(1)*(fisher(8)*(fisher(15)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))+fisher(33)*(fisher(16)*fisher(24)-fisher(22)*fisher(18)))-fisher(14)*(fisher(9)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))+fisher(33)*(fisher(10)*fisher(24)-fisher(22)*fisher(12)))+fisher(20)*(fisher(9)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))-fisher(15)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))+fisher(33)*(fisher(10)*fisher(18)-fisher(16)*fisher(12)))-fisher(32)*(fisher(9)*(fisher(16)*fisher(24)-fisher(22)*fisher(18))-fisher(15)*(fisher(10)*fisher(24)-fisher(22)*fisher(12))+fisher(21)*(fisher(10)*fisher(18)-fisher(16)*fisher(12))))-fisher(7)*(fisher(2)*(fisher(15)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))+fisher(33)*(fisher(16)*fisher(24)-fisher(22)*fisher(18)))-fisher(14)*(fisher(3)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(24)-fisher(22)*fisher(6)))+fisher(20)*(fisher(3)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))-fisher(15)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(18)-fisher(16)*fisher(6)))-fisher(32)*(fisher(3)*(fisher(16)*fisher(24)-fisher(22)*fisher(18))-fisher(15)*(fisher(4)*fisher(24)-fisher(22)*fisher(6))+fisher(21)*(fisher(4)*fisher(18)-fisher(16)*fisher(6))))+fisher(13)*(fisher(2)*(fisher(9)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))+fisher(33)*(fisher(10)*fisher(24)-fisher(22)*fisher(12)))-fisher(8)*(fisher(3)*(fisher(22)*fisher(36)-fisher(34)*fisher(24))-fisher(21)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(24)-fisher(22)*fisher(6)))+fisher(20)*(fisher(3)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))-fisher(9)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(12)-fisher(10)*fisher(6)))-fisher(32)*(fisher(3)*(fisher(10)*fisher(24)-fisher(22)*fisher(12))-fisher(9)*(fisher(4)*fisher(24)-fisher(22)*fisher(6))+fisher(21)*(fisher(4)*fisher(12)-fisher(10)*fisher(6))))-fisher(19)*(fisher(2)*(fisher(9)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))-fisher(15)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))+fisher(33)*(fisher(10)*fisher(18)-fisher(16)*fisher(12)))-fisher(8)*(fisher(3)*(fisher(16)*fisher(36)-fisher(34)*fisher(18))-fisher(15)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(18)-fisher(16)*fisher(6)))+fisher(14)*(fisher(3)*(fisher(10)*fisher(36)-fisher(34)*fisher(12))-fisher(9)*(fisher(4)*fisher(36)-fisher(34)*fisher(6))+fisher(33)*(fisher(4)*fisher(12)-fisher(10)*fisher(6)))-fisher(32)*(fisher(3)*(fisher(10)*fisher(18)-fisher(16)*fisher(12))-fisher(9)*(fisher(4)*fisher(18)-fisher(16)*fisher(6))+fisher(15)*(fisher(4)*fisher(12)-fisher(10)*fisher(6))))+fisher(31)*(fisher(2)*(fisher(9)*(fisher(16)*fisher(24)-fisher(22)*fisher(18))-fisher(15)*(fisher(10)*fisher(24)-fisher(22)*fisher(12))+fisher(21)*(fisher(10)*fisher(18)-fisher(16)*fisher(12)))-fisher(8)*(fisher(3)*(fisher(16)*fisher(24)-fisher(22)*fisher(18))-fisher(15)*(fisher(4)*fisher(24)-fisher(22)*fisher(6))+fisher(21)*(fisher(4)*fisher(18)-fisher(16)*fisher(6)))+fisher(14)*(fisher(3)*(fisher(10)*fisher(24)-fisher(22)*fisher(12))-fisher(9)*(fisher(4)*fisher(24)-fisher(22)*fisher(6))+fisher(21)*(fisher(4)*fisher(12)-fisher(10)*fisher(6)))-fisher(20)*(fisher(3)*(fisher(10)*fisher(18)-fisher(16)*fisher(12))-fisher(9)*(fisher(4)*fisher(18)-fisher(16)*fisher(6))+fisher(15)*(fisher(4)*fisher(12)-fisher(10)*fisher(6)))))/det_fish;
sigy_crlb = -(-(fisher(1)*(fisher(8)*(fisher(15)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))+fisher(27)*(fisher(16)*fisher(23)-fisher(22)*fisher(17)))-fisher(14)*(fisher(9)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))+fisher(27)*(fisher(10)*fisher(23)-fisher(22)*fisher(11)))+fisher(20)*(fisher(9)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))-fisher(15)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))+fisher(27)*(fisher(10)*fisher(17)-fisher(16)*fisher(11)))-fisher(26)*(fisher(9)*(fisher(16)*fisher(23)-fisher(22)*fisher(17))-fisher(15)*(fisher(10)*fisher(23)-fisher(22)*fisher(11))+fisher(21)*(fisher(10)*fisher(17)-fisher(16)*fisher(11))))-fisher(7)*(fisher(2)*(fisher(15)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))+fisher(27)*(fisher(16)*fisher(23)-fisher(22)*fisher(17)))-fisher(14)*(fisher(3)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(23)-fisher(22)*fisher(5)))+fisher(20)*(fisher(3)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))-fisher(15)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(17)-fisher(16)*fisher(5)))-fisher(26)*(fisher(3)*(fisher(16)*fisher(23)-fisher(22)*fisher(17))-fisher(15)*(fisher(4)*fisher(23)-fisher(22)*fisher(5))+fisher(21)*(fisher(4)*fisher(17)-fisher(16)*fisher(5))))+fisher(13)*(fisher(2)*(fisher(9)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))+fisher(27)*(fisher(10)*fisher(23)-fisher(22)*fisher(11)))-fisher(8)*(fisher(3)*(fisher(22)*fisher(29)-fisher(28)*fisher(23))-fisher(21)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(23)-fisher(22)*fisher(5)))+fisher(20)*(fisher(3)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))-fisher(9)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(11)-fisher(10)*fisher(5)))-fisher(26)*(fisher(3)*(fisher(10)*fisher(23)-fisher(22)*fisher(11))-fisher(9)*(fisher(4)*fisher(23)-fisher(22)*fisher(5))+fisher(21)*(fisher(4)*fisher(11)-fisher(10)*fisher(5))))-fisher(19)*(fisher(2)*(fisher(9)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))-fisher(15)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))+fisher(27)*(fisher(10)*fisher(17)-fisher(16)*fisher(11)))-fisher(8)*(fisher(3)*(fisher(16)*fisher(29)-fisher(28)*fisher(17))-fisher(15)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(17)-fisher(16)*fisher(5)))+fisher(14)*(fisher(3)*(fisher(10)*fisher(29)-fisher(28)*fisher(11))-fisher(9)*(fisher(4)*fisher(29)-fisher(28)*fisher(5))+fisher(27)*(fisher(4)*fisher(11)-fisher(10)*fisher(5)))-fisher(26)*(fisher(3)*(fisher(10)*fisher(17)-fisher(16)*fisher(11))-fisher(9)*(fisher(4)*fisher(17)-fisher(16)*fisher(5))+fisher(15)*(fisher(4)*fisher(11)-fisher(10)*fisher(5))))+fisher(25)*(fisher(2)*(fisher(9)*(fisher(16)*fisher(23)-fisher(22)*fisher(17))-fisher(15)*(fisher(10)*fisher(23)-fisher(22)*fisher(11))+fisher(21)*(fisher(10)*fisher(17)-fisher(16)*fisher(11)))-fisher(8)*(fisher(3)*(fisher(16)*fisher(23)-fisher(22)*fisher(17))-fisher(15)*(fisher(4)*fisher(23)-fisher(22)*fisher(5))+fisher(21)*(fisher(4)*fisher(17)-fisher(16)*fisher(5)))+fisher(14)*(fisher(3)*(fisher(10)*fisher(23)-fisher(22)*fisher(11))-fisher(9)*(fisher(4)*fisher(23)-fisher(22)*fisher(5))+fisher(21)*(fisher(4)*fisher(11)-fisher(10)*fisher(5)))-fisher(20)*(fisher(3)*(fisher(10)*fisher(17)-fisher(16)*fisher(11))-fisher(9)*(fisher(4)*fisher(17)-fisher(16)*fisher(5))+fisher(15)*(fisher(4)*fisher(11)-fisher(10)*fisher(5)))))/det_fish);

crlbs = [ xf_crlb, yf_crlb, N_crlb, sigx_crlb, sigy_crlb, off_crlb];
fits = [beta0(1), beta0(2), beta0(3),beta0(4),beta0(5),beta0(6)];
end