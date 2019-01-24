%% psf_localization
% This script is to localize all molecules associated with the PSF_cut
% output. Frames are expected to come in a cell format with all frames
% associated with one region in a single cell element. It will loop over
% all cell elements

clearvars; close all; clc;
% load('PSF_2.mat');
[fname, fpath] = uigetfile();
load([fpath,fname]);
load('z_calib.mat');
pixw = 4;
q = 0.128;
for i = 1:numel(A)
    iprod = A{i};
    [m,n,o] = size(iprod);
%     tic
%     thrsh =  mean(iprod(:));
%     toc
%     thrsh = 20/100*mean(max(max(iprod)));
%     ip = iprod(:,:,3) - iprod(:,:,1);
    for k = 1:o
        dps(:,:,k) = get_das_peaks(iprod(:,:,k),max(max(iprod(:,:,k))));
    end
    sum(dps(:))
    
    clear ip ipf i1
%     dps(:,:,3) = dps;
%     dps(:,:,1) = 0*dps(:,:,1);
%     dps(:,:,2) = 0*dps(:,:,1);
    % divide up the data
    [iloc, fnum, cents] = divide_up(iprod, pixw, dps);
    
    % remove duplicate data
    [ind] = find_dupes(cents,fnum);
%     iloc(:,:,ind) = [];
%     cents(ind,:) = [];
%     fnum(ind) = [];
    % Localize the Data
    % [xf_all,xf_crlb, yf_all,yf_crlb,sigx_all, sigx_crlb, sigy_all, sigy_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs_sigs(iloc, fnum, cents, angle);
    % zf_all = getdz(sigx_all,sigy_all)/q;
    % [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs(iloc, fnum, cents, angle);zf_all = zf_all/q;                        % This is to handle Z informtation uncomment once calibration is fixed
%     [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(iprod, fnum, cents*0+floor(m/2)+1, cal, pixw);
    [fits, crlbs, llv, framenumber] = slim_locs(iloc, fnum, cents, cal.ang);
    zf = getdz(fits(:,4),fits(:,5),cal.z_cal)/q;
    coords = [fits(:,1:2),zf];
    [ncoords] = astig_tilt(coords,cal);
% [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(iloc, fnum, cents, cal, pixw);
    try
        zf_crlb = zf_crlb/(q)^2; % this puts the CRLB in units of pix^2
        zf_all = zf_all/q; % this puts zf_all in units of pix
        
        %% Find all molecules that pass the initial localization requirements
        
        
        
        % Save the Analysis
        %  save([an_dir,'\', fname(1:end-4),'_dast.mat'], 'zf_all','sigx_all' ,'sigy_all','sigx_crlb','sigy_crlb','y','iloc','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','pixw','q','pix2pho');
        save([fname(1:end-4),'_num_',num2str(i),'_dast.mat'], 'cents','zf_all','zf_crlb','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','iters','pixw','q');
    catch lsterr
    end
end