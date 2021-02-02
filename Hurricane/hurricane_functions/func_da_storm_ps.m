function func_da_storm_ps(fname, q, pix2pho, pixw,thresh, angle, sv_im, mi1, choices,splits)

% Convert Variabls
% pix2pho = single(pix2pho);
q = double(q);
try
    
    comp_name = get_computer_name();
    load([comp_name , '\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat']);
    load([comp_name , '\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat'], 'split');
    load([comp_name , '\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat']);

    i1 = (readtiff(fname) - mi1);

    [m,n,o] = size(i1);
%     i1 = i1(:,:,1:2:o);
    iprod = gpu_rball(i1);
    if choices(4) == 1
        writetiff(iprod,[data_d,'\Rolling_Ball\',fname(1:end-4),'_rb.tif']);
    end

    
    clear i1 mi1;
    % Wavelet deniosing for identification
    iwaves = gpu_waves(iprod);
    se = strel('Disk',1);
    ifind = imerode(iwaves,se);
    if choices(1) == 1
        writetiff(ifind,[data_d,'\Waves\',fname(1:end-4),'_waves.tif']);
    end
    

    if choices(5) == 1 % User intended to use dual channel w/ both colors

            ifind = func_image_block(ifind,split,1); 
    elseif choices(5) == 2 % 2 = orange only channel intended
        ifind = func_image_red_block(ifind,split);
    elseif choices(5) == 3 % 3 = red only channel intended
        ifind = func_image_orange_block(ifind,split);
    end

     %Identification
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
    % Fitting
    switch choices(5)
        case 1 % 2 color orange first
            split = 171; % fill in this value based off data
            id = cents(:,1) < split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'red');
            id = logical(1-id);
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'orange', cdata);
        case 2 % Red only localization
            cdata = fit_channel_data(iloc, fnum, cents, cal, 'red');
        case 3 % orange only localization
            cdata = fit_channel_data(iloc, fnum, cents, cal, 'orange');
        case 4 % calibration localization
            split = 171; % fill in this value based off data
            id = cents(:,1) < split; % red localization
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'red');
            id = logical(1-id);
            cdata = fit_channel_data(iloc(:,:,id), fnum(id), cents(id,:), cal, 'orange', cdata);
    end
    
    % save results
        save(['Analysis\', fname(1:end-9),'_dast.mat'],  'cdata', 'pixw','q','cal');

    
catch lsterr
    disp(lsterr.message)
    disp(lsterr.stack(1))
    
    imagesc(iprod(:,:,1))
    title('Odd Frame')
    figure
    
    imagesc(iprod(:,:,101))
    title('Even Frame')
    waitforbuttonpress
    
end