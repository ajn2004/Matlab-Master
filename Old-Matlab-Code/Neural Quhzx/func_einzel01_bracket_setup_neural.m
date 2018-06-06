function [bkgn] = func_einzel01_setup_bkgn_only(data_dir,base_name,an_dir,q,n_start,n_bkgn, pix_to_pho)
% Original version from Sam Hess lab
% Modified by AJN
global xpix ypix wbox;

imagefile = strcat(data_dir,base_name);

if(exist('bkgn','var')==0)
    bkgn=bg_noise_calc01(imagefile,pix_to_pho,n_bkgn,base_name);
        
    answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),'. Yes to continue, No to restart.']));
    if ~strcmp(answer,'Yes')
        return;
    end
end



