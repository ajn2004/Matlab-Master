clear all


i1 = uint16(2*ones(100,100));
varian = ones(100,100);
sig = 3;
for i = 15:20:75
    for j = 15:20:75
        i1(i:i+5,j:j+5)=5;
    end
end

for i = 15:20:75
    for j = 15:20:75
        i1(i:i+10,j:j+10)=i1(i:i+10,j:j+10)+10;
    end
end


noise1 = double(imnoise(uint16(i1),'poisson'));
noisyim = double(i1)+noise1;
% % subplot(1,3,1);imshow(noisyim,[min(noisyim(:)), max(noisyim(:))]);
% % subplot(1,3,2);imshow(i1,[min(i1(:)), max(i1(:))]);
% % subplot(1,3,3);imshow(i1,[min(noise1(:)), max(noise1(:))]);
% % rball=3; %radius of rolling ball
% se = strel('square',rball); %structural element, i.e. rolling ball
% FWHM=1; %FWHM of gaussian smoothing in pixels
% rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
% kw=10; %kernal width of smoothing function
% [Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
% kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
% gs=exp(-2*kd.*kd/(rk*rk));
% gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
% 
% 
% gs=exp(-2*kd.*kd/(rk*rk));
% gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
% 
%     clims=[min(i1(:)), max(i1(:))];    
%     
%     i1_gs = uint16(conv2(noise1,gs,'same')); %smoothed original
%     bkg = double(imopen(i1_gs,se));
%     
%     iprod=double(i1)-bkg;
%     iprod=iprod.*(iprod>0); %set negative values to 0
iprod = func_scmos_var_smoother(noise1,sig,varian);
figure
subplot(1,3,1);imshow(iprod,[min(iprod(:)),max(iprod(:))]);
title('Corrected noise1');
subplot(1,3,2);imshow(noise1,[min(noise1(:)),max(noise1(:))]);
title('bkg');
subplot(1,3,3);imshow(i1,[min(i1(:)), max(i1(:))]);
title(['original sig=',num2str(sig)]);

