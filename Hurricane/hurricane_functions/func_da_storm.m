function func_da_storm(fname, data_d, an_dir, q, pixw,thresh, mi1, choices)

% Convert Variabls
% pix2pho = single(pix2pho);
q = double(q);
try
    

    comp_name = get_computer_name();
    load([comp_name , '\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat']);
    try
        load( 'z_calib.mat');
    catch
        load([comp_name , '\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat']);
    end
    load([comp_name , '\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat'],'split');

    i1 = (readtiff(fname) - mi1);

    [m,n,o] = size(i1);
%     i1 = i1(:,:,1:2:o);
    iprod = gpu_rball(i1);
    if choices(4) == 1
        writetiff(iprod,[data_d,'\Rolling_Ball\',fname(1:end-4),'_rb.tif']);
    end

    
    clear i1 mi1;
    iwaves = gpu_waves(iprod);
    se = strel('Disk',1);
    ifind = imerode(iwaves,se);
    if choices(1) == 1
        writetiff(ifind,[data_d,'\Waves\',fname(1:end-4),'_waves.tif']);
    end
    

%     if choices(5) == 1 % User intended to use dual channel w/ both colors
%         if odd_red_percentage < even_red_percentage % Even ratio larger than odd indicates red molecules are on even channels
%             ifind = func_image_block(ifind,split,1);
%         else
%             ifind = func_image_block(ifind,split,2); 
%         end
%     elseif choices(5) == 2 % 2 = orange only channel intended
%         ifind = func_image_red_block(ifind,split);
%     elseif choices(5) == 3 % 3 = red only channel intended
%         ifind = func_image_orange_block(ifind,split);
%     end

    %Peak identification
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
    
    
    cal.q = q;
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

    % Fitting
    switch choices(5)
        case 1 % 2 color orange first
            split = 181; % fill in this value based off data
            id = cents(:,1) < split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'red');
            id = logical(1-id);
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'orange', cdata);
        case 2 % Red only localization
            split = 181; % fill in this value based off data
            id = cents(:,1) < split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'red');
        case 3 % orange only localization
            split = 181; % fill in this value based off data
            id = cents(:,1) > split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'orange');
        case 4 % calibration localization
            split = 181; % fill in this value based off data
            id = cents(:,1) < split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'red');
            id = logical(1-id);
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'orange', cdata);
    end
    
    % save results
    save([an_dir,'\', fname(1:end-4),'_dast.mat'],  'cdata', 'pixw','q','cal');

    
catch lsterr
    disp(lsterr)
    
    imagesc(iprod(:,:,1))
    title('Odd Frame')
    figure
    
    imagesc(iprod(:,:,101))
    title('Even Frame')
    waitforbuttonpress
    
end

end