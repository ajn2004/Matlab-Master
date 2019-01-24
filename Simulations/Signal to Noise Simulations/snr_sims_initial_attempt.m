% % Make a bunch of noise
% % This simulation will take latex bead PSFs and add amounts of noise to the
% % image and attempt to localize the result to compare the signal to noise
% % verusus localization precision
scale = 100;
iprod = A{1}/scale;
% for i = 1:o
%     iprod(:,:,i) = iprod(:,:,i)./(max(max(iprod(:,:,i))));
% end


% iprod = A{1};
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
%     [ind] = find_dupes(cents,fnum);
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
totfit{1} = fits;
totcord{1} = ncoords;
totcrl{1} = crlbs;
meanS = mean(fits(:,3));
snrs = [60, 80, 100, 200, 300,400, 500];
offsets = meanS^2./snrs.^2;
for i = 1:numel(offsets)
    iprod = A{1}/scale;
    for j = 1:o % add noise to each image
        noisefrm = double(imnoise(uint16(round(offsets(i)*ones(m,n))),'poisson'));
        noisefrm1 = double(imnoise(uint16(round(offsets(i)*ones(m,n))),'poisson'));
        iprod(:,:,j) = iprod(:,:,j) + noisefrm;
        iprod(:,:,j) = iprod(:,:,j).*(iprod(:,:,j)>0);
    end
    [iloc, fnum, cents] = divide_up(iprod, pixw, dps);
    figure
    surf(iprod(:,:,1))
    title(['SNR ',num2str(snrs(i))])
    % remove duplicate data
%     [ind] = find_dupes(cents,fnum);
%     iloc(:,:,ind) = [];
%     cents(ind,:) = [];
%     fnum(ind) = [];
    % Localize the Data
    % [xf_all,xf_crlb, yf_all,yf_crlb,sigx_all, sigx_crlb, sigy_all, sigy_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs_sigs(iloc, fnum, cents, angle);
    % zf_all = getdz(sigx_all,sigy_all)/q;
    % [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, y, inloc, xin, yin] = da_locs(iloc, fnum, cents, angle);zf_all = zf_all/q;                        % This is to handle Z informtation uncomment once calibration is fixed
%     [xf_all,xf_crlb, yf_all,yf_crlb,zf_all, zf_crlb, N, N_crlb,off_all, off_crlb, framenum_all, llv, iters] = da_splines(iprod, fnum, cents*0+floor(m/2)+1, cal, pixw);
    [fits, crlbs, llv, framenumber] = slim_locs(iloc, fnum, cents, cal.ang);
%     try
    zf = getdz(fits(:,4),fits(:,5),cal.z_cal)/q;
        coords = [fits(:,1:2),zf];
        try
        [ncoords] = astig_tilt(coords,cal);
        catch
            ncoords = [];
        end
        totfit{i+1} = fits;
totcord{i+1} = ncoords;
totcrl{i+1} = crlbs;
end