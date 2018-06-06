% Last modified by AJN 3/21/16
% a function version of einzel to choose background levels in
% bracket_quhzx_neural_1.m
function [bkgn] = func_einzel_2color05_setup(data_dir,base_name,n_bkgn,pix_to_pho)



imagefile = strcat(data_dir,base_name);

%background niose estimation
if(exist('bkgn','var')==0)
    bkgn1=bg_noise_calc01(imagefile,pix_to_pho,n_bkgn, base_name);
    bkgn2=bg_noise_calc01(imagefile,pix_to_pho,n_bkgn, base_name);
    bkgn=sqrt(bkgn1^2+bkgn2^2);
    
    answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),'. Yes to continue, No to restart.']));
    if ~strcmp(answer,'Yes')
        return;
    end
end

