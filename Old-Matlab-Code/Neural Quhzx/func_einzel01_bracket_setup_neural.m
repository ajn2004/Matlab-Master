function [bkgn] = func_einzel01_setup_bkgn_only(data_dir,base_name,an_dir,q,n_start,n_bkgn, pix_to_pho)

global xpix ypix wbox;

imagefile = strcat(data_dir,base_name);
% rbox=3;
% 
% wvlnth=580/1000; %convert wavelength from nm to um       
% NA=1.2;
% psf_scale=1.2;
% box_overlap_factor = 1.5; %if center to center closer, don't include either
% 
% w_mask = round(rbox*box_overlap_factor); %width of masking when a "high" pixel is identified
% wbox=2*rbox+1;
% [xpix0,ypix0] = meshgrid(-2*rbox:2*rbox,-2*rbox:2*rbox);
% [xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
% 
% psf_w0 = psf_scale*0.55*wvlnth/NA/1.17; 
% psf_std=psf_w0/2; %standard deviation of psf
% psf_w02=(psf_w0/q)*(psf_w0/q); %square of 1/e^2 radius in pixels
% 
% rball=6; %radius of rolling ball
% se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball
% FWHM=1; %FWHM of gaussian smoothing in pixels
% rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
% kw=20; %kernal width of smoothing function
% [Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
% kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
% gs=exp(-2*kd.*kd/(rk*rk));
% gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
% 
% %initialize arrays
% n_init=500000;
% xcm_all=zeros(n_init,1);
% ycm_all=zeros(n_init,1);
% xf_all=zeros(n_init,1);
% yf_all=zeros(n_init,1);
% a0_all=zeros(n_init,1);
% r0_all=zeros(n_init,1);
% off_all=zeros(n_init,1);
% framenum_all=zeros(n_init,1);
% xf_err_all=zeros(n_init,1);
% yf_err_all=zeros(n_init,1);
% a0_err_all=zeros(n_init,1);
% r0_err_all=zeros(n_init,1);
% off_err_all=zeros(n_init,1);
% grab_sum_all=zeros(n_init,1);
% 
% total_molecules=0;
% n_fail_a0=0;
% n_fail_outbox=0;



% bkgn=3.1168;% background noise in photons, comment this out to auto calculate
%auto calc bkgnd noise if bkgn is commented out above
if(exist('bkgn','var')==0)
    bkgn=bg_noise_calc01(imagefile,pix_to_pho,n_bkgn,base_name);
        
    answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),'. Yes to continue, No to restart.']));
    if ~strcmp(answer,'Yes')
        return;
    end
end



% fileloop = n_start;
% while fileloop>0
% 
% 
%     drawnow;
%     
%     i1=double(imread(imagefile,fileloop))/pix_to_pho;
%     
%     [yw,xw]=size(i1);
%         
%     i1_gs = uint16(conv2(i1,gs,'same')); %smoothed original
%     bkg = double(imopen(i1_gs,se));
%     
%     iprod=i1-bkg;
%     iprod=iprod.*(iprod>0); %set negative values to 0
%     
%     %%%%%%%%%
% %     bg_RBSub_noise=bg_noise_image(iprod);   
% %     answer = questdlg(sprintf(['RB subtracted image Noise: ',num2str(bg_RBSub_noise,'%2.2f'),', Thresh: ',num2str(min_thresh),'. Yes to continue, No to restart.']));
% %     if ~strcmp(answer,'Yes')
% %         return;
% %     end
%     %%%%%%%%%%%%
%     
%     high_pixel_mask = zeros(yw,xw);
%     n_boxes=0;
%     boxes_xy=zeros(10000,3);
% 
%     [iy,ix] = find(iprod >= min_thresh); %find all pixels in iprod above the threshold
% 
%     high_pix_inframe = size(iy,1); %the number of pixels above thresh in frame (iy:ix)
% 
%     while true
%         pix_val = 0;
%         for y1=1:1:high_pix_inframe
%             for x1=1:1:high_pix_inframe
%                 if (iprod(iy(y1),ix(x1))>pix_val && high_pixel_mask(iy(y1),ix(x1))==0 ...
%                         && iy(y1)<yw-rbox-1 && iy(y1)>rbox+1 && ix(x1)<xw-rbox-1 && ix(x1)>rbox+1)
%                    pix_val = iprod(iy(y1),ix(x1));  
%                    high_pixel_y = iy(y1);
%                    high_pixel_x = ix(x1);
%                 end
%             end
%         end
%         if pix_val < min_thresh %break the loop if no new high pixels were found
%             break
%         end
% 
%         x0_box=high_pixel_x-rbox;
%         y0_box=high_pixel_y-rbox;
%         x1_box=high_pixel_x+rbox;
%         y1_box=high_pixel_y+rbox;
% 
%         x0_mask=high_pixel_x-w_mask;
%         if x0_mask < 1
%             x0_mask = 1;
%         end
%         x1_mask=high_pixel_x+w_mask;
%         if x1_mask > xw
%             x1_mask = xw;
%         end
%         y0_mask=high_pixel_y-w_mask;
%         if y0_mask < 1
%             y0_mask = 1;
%         end
%         y1_mask=high_pixel_y+w_mask;
%         if y1_mask > yw
%             y1_mask = yw;
%         end
% 
%         high_pixel_mask(y0_mask:y1_mask,x0_mask:x1_mask)=1;
% 
%         grab=iprod(y0_box:y1_box,x0_box:x1_box);
%         grab_sum=sum(sum(grab));
% 
%         %calculate C.O.M. for intial guess in fit
%         xm_sum=0;
%         ym_sum=0;
%         m_sum=0;
% 
%         for i=x0_box:x1_box
%             for j=y0_box:y1_box
%                 xind=floor(i);
%                 yind=floor(j);
%                 intens=iprod(yind,xind);
%                 xm_sum=xm_sum+xind*intens;
%                 ym_sum=ym_sum+yind*intens;
%                 m_sum=m_sum+intens;
%             end
%         end
% 
%         x_cm=xm_sum/m_sum;
%         y_cm=ym_sum/m_sum;
% 
%         xc_box=(x0_box+x1_box)*0.5;
%         yc_box=(y0_box+y1_box)*0.5;
% 
%         xguess=x_cm-xc_box;
%         yguess=y_cm-yc_box;
% 
%         for i=1:wbox
%             for j=1:wbox
%                 k=(i-1)*wbox+j;
%                 xymerge(k)=0;
%                 zmerge(k)=grab(i,j);
%             end
%         end
% 
%         beta0=[xguess,yguess,50,psf_w0/q,min(grab(:))]; % x,y,a0,r0,offset
%         [betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);
%         ci = nlparci(betafit,resid,'covar',COVB); %calculate error estimates on parameters
%         ci_err=(ci(:,2)-ci(:,1))/2;
% 
%         yf=betafit(2)+yc_box;
%         xf=betafit(1)+xc_box;
%         a0=betafit(3);
%         r0=abs(betafit(4));
%         off=betafit(5);
%         
%         failed=0; %flag for failed localization
%         if(a0 < 0)
%             n_fail_a0=n_fail_a0+1;
%             failed=1;
%         end
%         if(xf > x1_box || xf < x0_box || yf > y1_box || yf < y0_box)
%             n_fail_outbox=n_fail_outbox+1;
%             failed=1;
%         end
% 
%         if(failed==0) %assign if fit criteria is satisfied
%             total_molecules=total_molecules+1;
%             xcm_all(total_molecules)=x_cm;
%             ycm_all(total_molecules)=y_cm;
%             xf_all(total_molecules)=xf;
%             yf_all(total_molecules)=yf;
%             a0_all(total_molecules)=a0;
%             r0_all(total_molecules)=r0;
%             off_all(total_molecules)=off;
%             framenum_all(total_molecules)=fileloop;
%             xf_err_all(total_molecules)=ci_err(1);
%             yf_err_all(total_molecules)=ci_err(2);
%             a0_err_all(total_molecules)=ci_err(3);
%             r0_err_all(total_molecules)=ci_err(4);
%             off_err_all(total_molecules)=ci_err(5);
%             grab_sum_all(total_molecules)=grab_sum;
%                        
%             n_boxes=n_boxes+1;
%             boxes_xy(n_boxes,1)=high_pixel_x;
%             boxes_xy(n_boxes,2)=high_pixel_y;
%             boxes_xy(n_boxes,3)=1;
%         end
%     end
%         clims = [min(iprod(:)), max(iprod(:))+2];
%         figure;
%         imagesc(iprod,clims), axis image, colormap gray, colorbar
%         title(['Frame: ' num2str(fileloop) ', ' num2str(total_molecules) ' molecules']);
%         hold on
%         draw_boxes(n_boxes,boxes_xy,rbox);
%         hold off
%         drawnow;
%         disp(['Threshold is currently ',num2str(min_thresh)]);
%         thresh_ok = input('Are you happy with the threshold value? (Y)es or (N)o or (S)kip','s');
%         if strcmp(thresh_ok,'Y')||strcmp(thresh_ok,'y') 
%             anlyze = input('Analyze image set using this value?(no to see another frame) ','s');
%             if strcmp(anlyze,'y')||strcmp(anlyze,'Y')
%                 real_thresh =min_thresh;
%                 fileloop= -1;
%                 close all
%             else
%                 fileloop=fileloop+1;
%                 close all
%             end   
%         elseif strcmp(thresh_ok,'S')||strcmp(thresh_ok,'s')
%             fileloop=fileloop+1;
%             close all
%         else
%              min_thresh = input('New Threshold Value? ');
%              close all
%         
%         end
%         
%     end
% end
% 
% 
