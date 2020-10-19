function func_da_storm(fname, data_d, an_dir, q, pixw,thresh, mi1, choices)

% Convert Variabls
% pix2pho = single(pix2pho);
q = double(q);
try
    
    try
    load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
    load('C:\Users\andre\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat')
    catch
        load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
        load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Hurricane\hurricane_functions\z_calib.mat')
    end

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
    
%     % Beginning of comment block for switching behavior
% %     % automatically detect switcher behavior
    dps = cpu_peaks(ifind(:,:,1:100),50,pixw);
    [iloc, fnum, cents] = divide_up(iprod(:,:,1:100), pixw, dps);
    
    % The ongoing pursuit of how to automatically detect which frame flips
    % first. This should be corrected in arduino code.
    if isempty(cents)
            dps = cpu_peaks(ifind(:,:,1:100),20,pixw);
            [iloc, fnum, cents] = divide_up(iprod(:,:,1:100), pixw, dps);
    end
    odd_red_percentage = sum(mod(fnum,2) == 1 & cents(:,1) < 180)/(sum(mod(fnum,2) == 1 & cents(:,1) > 180)+1); % on odd frames, ratio of red / orange molecules
    even_red_percentage = sum(mod(fnum,2) == 0 & cents(:,1) < 180)/(sum(mod(fnum,2) == 0 & cents(:,1) > 180)+1); % on even frames, ratio of red / orange molecules
    
%     % Count percentage of molecules in even and odd channels
%     ind = mod(fnum,2) == 0;
%     red_evens = sum(cents(ind,1)<180)/sum(ind);
%     orange_evens = sum(cents(ind,1)>180)/sum(ind);
%     even_red_brights = sum(sum(iloc(:,:,cents(ind,1)<180)));
%     even_orange_brights = sum(sum(iloc(:,:,cents(ind,1)>180)));
%     even_ratio = mean(even_red_brights)/mean(even_orange_brights); % shows average 'brightness' of red to orange
%     
%     ind = mod(fnum,2) == 1;
%     red_odds = sum(cents(ind,1)<180)/sum(ind);
%     orange_odds = sum(cents(ind,1)>180)/sum(ind);
%     odd_red_brights = sum(sum(iloc(:,:,cents(ind,1)<180)));
%     odd_orange_brights = sum(sum(iloc(:,:,cents(ind,1)>180)));
%     odd_ratio = mean(odd_red_brights)/mean(odd_orange_brights);

    if choices(5) == 1 % User intended to use dual channel w/ both colors
        if odd_red_percentage < even_red_percentage % Even ratio larger than odd indicates red molecules are on even channels
            ifind = func_image_block(ifind,split,1);
        else
            ifind = func_image_block(ifind,split,2); 
        end
    elseif choices(5) == 2 % 2 = orange only channel intended
        ifind = func_image_red_block(ifind,split);
    elseif choices(5) == 3 % 3 = red only channel intended
        ifind = func_image_orange_block(ifind,split);
    end
% % End of switching behavior comment block

%     % If we're doing 2 color, block out frame we're not interested in
%     if choices(5) == 1 % one equl dual color switcher
%         load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
%         ifind = func_image_block(ifind,split);
%     elseif choices(5) == 2 % 2 = orange only channel
%         load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
%         ifind = func_image_red_block(ifind,split);
%     elseif choices(5) == 3 % 0 = red only channel
%         load('C:\Users\andre\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
%         ifind = func_image_orange_block(ifind,split);
%     end
   
    
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
        
        zf = getdz(abs(fits(:,4)),abs(fits(:,5)),cal.red.z_cal,2)/q;
        coords = [fits(:,1:2),zf];
        [ncoords] = astig_tilt(coords,cal.red);
        save([an_dir,'\', fname(1:end-4),'_dast.mat'],  'pixw','q','ncoords','fits','crlbs','llv','framenumber','cal');
    end
    id = cents(:,1) < split; % Identify localizations below the split
    if choices(5) == 1 || choices(5) ==3
        %     load('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\2-Channel Codes\2_color_calibration.mat', 'split', 'o2rx','o2ry');
        
        %% First fit is all red, so those can be immediately
        if sum(id)>0
            
            [fits, crlbs, llv, framenumber] = slim_locs(iloc(:,:,id), fnum(id), cents(id,:), cal.red.ang);
            
            % As everywhere in the equations used sigma is squared, we can without
            % loss of generality make these fits positive definite
            fits(:,4) = abs(fits(:,4));
            fits(:,5) = abs(fits(:,5));
            
            % Z calculations
            %             zf = getdz(cdata.red.fits(:,4),cdata.red.fits(:,5),cal.red.z_cal,2)/q; % Z assignment, this variable will be updated and stored elsewhere
            %             ncoords = astig_tilt([cdata.red.fits(:,1:2),zf],cal.red); % corrections due to astigmatism
            zf = get_spline_z(fits(:,4),fits(:,5),cal.red); % New z_registration based off spline 3d calibration
            % remove failed z assignments
            ind = zf == -5;
            fits(ind,:) = [];
            crlbs(ind,:) = [];
            llv(ind) = [];
            framenumber(ind) = [];
            zf(ind) = [];
            
            
            % Put data into cdata structure
            for i = 1:6
                cdata.red.fits(:,i) = fits(:,i);
                cdata.red.crlbs(:,i) = crlbs(:,i);
            end
            cdata.red.llv = llv;
            cdata.red.framenumber = framenumber;
            
            ncoords = make_astigmatism_corrections([cdata.red.fits(:,1:2),zf/q],cal.red,q);
            % Assign fixed coordinates which are in microns at this point
            cdata.red.xf = ncoords(:,1);
            cdata.red.yf = ncoords(:,2);
            cdata.red.zf = ncoords(:,3);
            
            red_index = cdata.red.xf  < 0 | cdata.red.yf <0 ;
            field_names = fieldnames(cdata.red);
            
            if sum(red_index)>0
                cdata.red.fits(red_index,:) = [];
                cdata.red.crlbs(red_index,:) = [];
            end
            
            for k=3:numel(field_names)
                if sum(red_index) > 0
                    cdata.red.(field_names{k})(red_index) = [];
                end
                
            end
            
            clear fits crlbs llv framenumber
        end
        %% Repeat above for orange
        id = logical(1-id); % Changes 0 -> 1 and 1 -> 0 flipping the ID so now we can fit orange
        if sum(id) >0 && (choices(5) == 1 || choices(5) == 2)
            
            [fits, crlbs, llv, framenumber] = slim_locs(iloc(:,:,id), fnum(id), cents(id,:), cal.orange.ang);
            % As everywhere in the equations used sigma is squared, we can without
            % loss of generality make these fits positive definite
            fits(:,4) = abs(fits(:,4));
            fits(:,5) = abs(fits(:,5));
            % Remove failed Z fits
            zf = get_spline_z(fits(:,4),fits(:,5),cal.orange); % New z_registration based off spline 3d calibration            
            ind = zf == -5;
%             fits(ind,:) = [];
%             crlbs(ind,:) = [];
%             llv(ind) = [];
%             framenumber(ind) = [];
%             zf(ind) = [];
            % Put data into cdata structure
            for i = 1:6
                cdata.orange.fits(:,i) = fits(:,i);
                cdata.orange.crlbs(:,i) = crlbs(:,i);
            end
            cdata.orange.llv = llv;
            cdata.orange.framenumber = framenumber;
            
            % Z calculations
            %zf = getdz(cdata.orange.fits(:,4),cdata.orange.fits(:,5),cal.orange.z_cal,2)/q; % Z assignment, this variable will be updated and stored elsewhere
            
            %ncoords = astig_tilt([cdata.orange.fits(:,1:2),zf],cal.orange); % corrections due to astigmatism
            %ncoords = astig_tilt([cdata.orange.fits(:,1:2),zf],cal.orange); % corrections due to astigmatism
            ncoords = make_astigmatism_corrections([cdata.orange.fits(:,1:2),zf/q],cal.orange,q);
            vec = xy_feature(ncoords(:,1),ncoords(:,2));
%             x = o2rx.'*vec.';
%             y = o2ry.'*vec.';
            % Assign fixed coordinates ncoords are in microns
            cdata.orange.xf = ncoords(:,1);
            cdata.orange.yf = ncoords(:,2);
            cdata.orange.zf = ncoords(:,3);
            % ensure transform is assigned to calibration variable
            cal.o2rx = o2rx;
            cal.o2ry = o2ry;
            field_names = fieldnames(cdata.orange);
        
            orange_index = cdata.orange.xf  < 0 | cdata.orange.yf <0 ;
            if sum(orange_index)>0
                cdata.orange.fits(orange_index,:) = [];
                cdata.orange.crlbs(orange_index,:) = [];
            end
            for k=3:numel(field_names)
                if sum(orange_index) > 0
                    cdata.orange.(field_names{k})(orange_index) = [];
                end
            end
        end
        % Remove problem entries
        

        
        
        
          
        
        save([an_dir,'\', fname(1:end-4),'_dast.mat'],  'cdata', 'pixw','q','cal');
    end
    
catch lsterr
    disp(lsterr)
    
    imagesc(iprod(:,:,1))
    title('Odd Frame')
    figure
    
    imagesc(iprod(:,:,101))
    title('Even Frame')
    waitforbuttonpress
    
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


% end
% catch lsterr
%      save([an_dir,'\', fname(1:end-4),'_dast.mat'], 'zf_all','sigx_all' ,'sigy_all','sigx_crlb','sigy_crlb','y','iloc','xf_all' , 'xf_crlb' , 'yf_all' , 'yf_crlb' , 'N' , 'N_crlb' ,'off_all' , 'off_crlb', 'framenum_all', 'llv','pixw','q','pix2pho');
end