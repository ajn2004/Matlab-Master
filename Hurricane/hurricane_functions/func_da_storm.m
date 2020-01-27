function func_da_storm(fname,data_d, an_dir, q, pix2pho, pixw,thresh, angle, sv_im, mi1, choices)

% Convert Variabls
% pix2pho = single(pix2pho);
q = double(q);
load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat')
% if exist([data_d, 'z_calib.mat'])
%     cal = load([data_d 'z_calib.mat']);
% else
%     cal = load('z_calib.mat');
% end
% mi1 = 0
% Load file and dark current background subtraction
i1 = (readtiff(fname) - mi1);
% i1 = sum(i1,3);
% i1 = i1.*(i1>0);
[m,n,o] = size(i1);
% i1(1:30,:,:) = 0;
% i1(m-30:m,:,:) = 0;
% i1(:,1:30,:) = 0;
% i1(:,n-30:n,:) = 0;
% Rolling Ball Background Subtract
% iprod = rollingball(i1);
iprod = gpu_rball(i1);
if choices(4) == 1
    writetiff(iprod,[data_d,'\Rolling_Ball\',fname(1:end-4),'_rb.tif']);
end
% iprod = bp_subtract(i1);
% iprod = imgaussfilt(i1,0.8947);
% iprod = i1;
% Peak Detection


% thrsh = 300/pix2pho;
% diprod = diff(iprod,3);
% for i = 1:o
% ifind = denoise_psf(iprod,2);
iwaves = gpu_waves(iprod);
se = strel('Disk',1);
ifind = imerode(iwaves,se);
if choices(1) == 1
    writetiff(ifind,[data_d,'\Waves\',fname(1:end-4),'_waves.tif']);
end
% thrsh = thresh/100*mean(max(max(iprod)));
% tic
% thrsh = 3*std(iprod(:)) + mean(iprod(:));
% toc
% thrsh = thresh/100*mean(max(max(diprod)));
% surf(max(ifind,[],3));
% thrsh = input('What should the threshold be? ');
% thrsh = min(iprod(:))*thresh/100;
% excerpt = 1;
dps = cpu_peaks(ifind,thresh,pixw);
if choices(2) == 1
    in_d_eye(iprod, dps, pixw);
end

clear ip ipf i1

% divide up the data

[iloc, fnum, cents] = divide_up(iprod, pixw, dps);
[m,n,o] = size(iloc);
% remove duplicate data
[ind] = find_fm_dupes(cents,fnum,pixw*1.5);
iloc(:,:,ind) = [];
cents(ind,:) = [];
fnum(ind) = [];



% Localize the Data
% [xf_all,xf_crlb, yf_all,yf_crlb,sigx_all, sigx_crlb, sigy_all, sigy_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs_sigs(iloc, fnum, cents, angle);
% zf_all = getdz(sigx_all,sigy_all)/q;
% [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs(iloc, fnum, cents, angle);zf_all = zf_all/q;                        % This is to handle Z informtation uncomment once calibration is fixed
% [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters, cex, cey] = da_splines(iloc, fnum, cents, cal, pixw);
% [~,~, ~,~,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters, cex, cey] = da_splines(iloc, fnum, cents, cal, pixw);
% i2 = reshape(iloc,m*n,o);
% save('thisbit.mat','iloc','cents','fnum','cal');
if choices(3) == 1
    writetiff(iloc,[data_d,'\psfs\',fname(1:end-4),'_psfs.tif']);
end
if choices(5) == 0
[fits, crlbs, llv, framenumber] = slim_locs(iloc, fnum, cents, cal.red.ang);
fits(:,4) = abs(fits(:,4));
fits(:,5) = abs(fits(:,5));

    zf = getdz(abs(fits(:,4)),abs(fits(:,5)),cal.z_cal)/q;
    coords = [fits(:,1:2),zf];
    [ncoords] = astig_tilt(coords,cal);
else
    load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat');
    id = cents(:,1) < split; % Identify localizations below the split
    % First fit is all red, so those can be immediately 
    [fits, crlbs, llv, framenumber] = slim_locs(iloc(:,id), fnum(id), cents(id,:), cal.red.ang);
    
    fits(:,4) = abs(fits(:,4));
    fits(:,5) = abs(fits(:,5));
    cdata = channel_correct(fits(:,1:2)); %% Identify and correct channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% and below
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% needs fixing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now
    % The idea is going to be to fit the 2 channels with different
    % calibration files for 3D to give us the best fit possible.
    try % Attempt 3D with 2 color calibration, if it doesn't work use 1 color calibration for both channels
        cdata.red.zf = getdz(abs(fits((cdata.id == 1),4)),abs(fits((cdata.id == 1),5)),cal.red.z_cal)/q;
        cdata.orange.zf = getdz(abs(fits((cdata.id == 2),4)),abs(fits((cdata.id == 2),5)),cal.orange.z_cal)/q;
    catch % this shouldn't matter after proper debugging
        cdata.red.zf = getdz(abs(fits((cdata.id == 1),4)),abs(fits((cdata.id == 1),5)),cal.z_cal)/q;
        cdata.orange.zf = getdz(abs(fits((cdata.id == 2),4)),abs(fits((cdata.id == 2),5)),cal.z_cal)/q;
    end
    % Organize fitting data with molecule type
    cdata.red.N = fits(cdata.id == 1,3);
    cdata.red.O = fits(cdata.id == 1,6);
    cdata.red.sx = fits(cdata.id == 1,4);
    cdata.red.sy = fits(cdata.id == 1,5);
    cdata.red.crlbs = crlbs(cdata.id == 1);
    cdata.red.llv = llv(cdata.id == 1);
    
    cdata.orange.N = fits(cdata.id == 2,3);
    cdata.orange.O = fits(cdata.id == 3,6);
    cdata.orange.sx = fits(cdata.id == 2,4);
    cdata.orange.sy = fits(cdata.id == 2,5);
    cdata.orange.crlbs = crlbs(cdata.id == 2);
    cdata.orange.llv = llv(cdata.id == 2);
end


% save('results_of_bump.mat','fnum','q','iloc','cal','cents');

% save('for_trial.mat','iloc'

% Save the Analysis
%  save([an_dir,'\', fname(1:end-4),'_dast.mat'], 'zf_all','sigx_all' ,'sigy_all','sigx_crlb','sigy_crlb','y','iloc','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','pixw','q','pix2pho');
% if strcmp(sv_im,'Y') || strcmp(sv_im,'y')
% save([an_dir,'\', fname(1:end-4),'_dast.mat'], 'cents','zf_all','zf_crlb','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','iters','pixw','q','pix2pho','ilocs');
% else
% figure
% imagesc(mean(iprod,3));
% hold on
% plot(fits(:,1),fits(:,2),'rx')
% hold off
% colormap('gray');

save([an_dir,'\', fname(1:end-4),'_dast.mat'],  'cdata', 'pixw','q','ncoords','fits','crlbs','llv','framenumber','cal');
% end
% catch lsterr
%      save([an_dir,'\', fname(1:end-4),'_dast.mat'], 'zf_all','sigx_all' ,'sigy_all','sigx_crlb','sigy_crlb','y','iloc','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','pixw','q','pix2pho');
end

