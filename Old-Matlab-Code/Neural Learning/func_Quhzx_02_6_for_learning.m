function [savename] = func_Quhzx_02_6_for_learning(fname, an_dir, iprod,bkgn,Theta1, Theta2, sigma2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func_Quhzx v 2.04
%
% Function compatibility is set up to be used with a bracketing program to
% go through all files in a data set
%
% Quhzx is a MLE based localization code capable of localizing molecular
% images from an emccd camera by fitting 2 dimensional Gaussian functions
% to an input image.
%
% The input of Quhzx is a .tif image or image stack.
% The output of Quhzx is a '.mat' file saved to a folder of the user's
% choice.
%
% Comments should help explain the function of Quhzx throughout the code
%
% This code was built based off the supplemental information in
% Smith et al. Nat Meth 2010
%
% This code was also heavily inspired by the einzel.m codes written by the
% Sam Hess Lab
%
% AJN 4/11/15  email: andrew.j.nelson@maine.edu
%
%
% v 1.02 replaced peak finding algorithm with more sensistive and faster
% algorithm
% Added Cramer Lao lower bound calculations
% v 1.03 preGPU release
% made major adjustments to flow of data in the code to be compatible with
% the mex gpu protocol.
% v 2.00  Quhzx has been crudely interfaced with the GPU
% v 2.01 Interfacing with GPU has been improved to handle variable image
% size without crashing. The log likelihood value and crlbs have been
% verified and can be relied on
% v 2.02 Attempting to minimize image loading time. By preallocating the
% variable i1 and passing the imfinfo variable to imread we are able to get
% ~17x speed increase in the loading of the file
% v 2.03 integrated with neural network
% V 2.04 integrated with gpu based neural network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User controlled variables
% min_thresh = 30;                 % minimum threshold for identifying a molecule. The units are photons this value will be updated during the code
% q = 0.152;                      % this is the pixel size in microns
% wvlngth = 573;                  % peak wavelength of molecular emission
% an_dir = 0;                   % analysis directory to where '.mat' files will be saved set to 0 to save in data directory
% n_start = 1;                   % starting frame
% n_end = 200;                      % ending frame type 0 to go to the last frame in your folder
% n_bkgn = 10;                    % frame to measure background
% pix2pho = 40;                   % pixel to photon ratio
% NA = 1.45;                      % Numerical aperture of the objective you used
% bkgn = 0;                       % background noise if known set to 0 to measure



%% File acquisition and loading
% [fname, fpath] = uigetfile('*.tif','Please select an tiff stack');          % user selects file to be analyzed
% [neur, neur_path] = uigetfile('*.mat','Please select a neural parameter file');
% load([neur_path,neur]);
% imagfo = imfinfo([fpath,fname]);                                          % variable to measure image info
% current_dir = pwd;
% cd(fpath);
% addpath(current_dir);
if an_dir ==0
    an_dir = fpath;
end
%% Measure background
% if bkgn == 0
%     check1 = 0;
%     while check1 == 0                                                      % while loop should allow for repeated measures until the user is satisfied
%         bkgn = bg_noise_calc02([fpath,fname],pix2pho,n_bkgn);
%
%         button = questdlg(['The background noise was measured to be ', num2str(bkgn),' photons. Use this value? (No to remeasure)']);  % really cool function
%         if strcmp(button,'Yes')
%             check1 = 1;
%         end
%     end
% end
% tic


%% setup localization variables

[xgrid, ygrid] = meshgrid(-ceil(sigma2)-1:ceil(sigma2)+1,-ceil(sigma2)-1:ceil(sigma2)+1);   % grids used for building fitting gaussians
wgrid = 2*ceil(sigma2)+3;                                                                   % width of the grids for fitting

i6 = [];
xf_all = [];
yf_all = [];
N = [];
off_all = [];
xf_crlb = [];
yf_crlb = [];
N_crlb = [];
off_crlb = [];
framenum_all = [];
llv = [];
i5 = [];
xc = [];
yc = [];
sigx_all =[];
sigy_all = [];
sigx_crlb = [];
sigy_crlb = [];
i6 = [];

runs = 1;

threadcount = 250;

% calculate memory requirements
rows = numel(iprod(:,1,1));
cols = numel(iprod(1,:,1));
ims = numel(iprod(1,1,:));

mem_size = rows*cols*ims*8;

max_neur = 200*200*500*8; %this number is based off simulations
a1 =[];
num_chunk = ceil(mem_size/max_neur);
chunk_size = ceil(ims/num_chunk);
neur_wait = waitbar(0, 'Neural Net calcualtion is 0% complete');
for i = 1:num_chunk
    if i~= num_chunk
        [atemp] = image_neural_2(iprod(:,:,1+ (i-1)*chunk_size: i*chunk_size), Theta1.', Theta2.', chunk_size);
        if i ==1
            a1 = atemp;
            clear atemp
        else
            a1 = cat(3,a1,atemp);
            clear atemp
        end
    else
        [atemp] = image_neural_3(iprod(:,:,1+ (i-1)*chunk_size:end), Theta1.', Theta2.', numel(iprod(1,1,1+ (i-1)*chunk_size:end)));
        a1 = cat(3,a1,atemp);
        clear atemp
    end
    waitbar(i/num_chunk,neur_wait,['Neural Net calcualtion is ', num2str(100*i/num_chunk),'% complete']);
end
disp('Neural calc done');
clear i1
close(neur_wait)
% close(neur_wait)
%% Frame Analysis Looping
disp('Localizing, please wait');
mol_num = 0;                                    % variable to keep track of molecule number
runs = 1;
% figure('units','Normalized','Outerposition',[0 0 1 1])
look = randi(numel(iprod(1,1,:)));
n_start = 1;
n_end = numel(iprod(1,1,:));
for framenum = 1:numel(iprod(1,1,:))                               % start on frame n_start and cycle through frame n_end
    
    clear i2 i3 cent xpeak ypeak num_high_pix removes
    %     waitbar(framenum/(n_end-n_start+1),hwait,[num2str(100*framenum/(1+n_end-n_start)),'% complete']);     % pdate waitbar
    i2 = iprod(:,:,i);   % frame to be analyzed is loaded
    
    %     [xpeak, ypeak] = grab_peaks(i2, min_thresh);                      % find all local maxima in the image
    [mol_locs, num_mol] = neural_gpu(a1(:,:,i));   % Perform neural net calculation
    if 1 - isempty(mol_locs)
        xpeak = mol_locs(:,1);
        ypeak = mol_locs(:,2);
    else
        xpeak = [];
        ypeak = [];
    end
    num_high_pix = numel(xpeak);                   % number of high pixels
    [yw, xw] = size(i2);                            % determines width of image
    
    
    %% Image Segmentation and localization identification
    % at this point only peaks above the minimum threshold will be recorded
    if 1 - isempty(xpeak)                  % thisif statement will skip blank frames
        if framenum == look
%             subplot(1,2,1);
            imagesc(i2);
            colormap('gray')
            axis image
            boxes([xpeak,ypeak],3.5,'g')
            title('Green boxes show identified regions')
            drawnow
            break
        end
        num_high_pix = numel(xpeak);
        m=1;                                    % initialize counting variable
        removes = 0;                        % initialize remove index
        for j = 1:num_high_pix                % cycles through all molecules except the very last
            if ~ismember(j,removes) && j < num_high_pix            % as long as current molecule hasn't already been removed
                for k = j+1:num_high_pix        % cycles through all other molecules
                    closeflag = 0;              % alarm variable is reset to 0
                    distance = ((xpeak(j)-xpeak(k))^2 + (ypeak(j)-ypeak(k))^2)^0.5;    % check distance of molecule pair
                    if distance < sigma2                                    % trips on pair that are closer than the diffraction limit
                        closeflag = 1;          % alarm variable triggered
                        break                   % stop the search
                    end
                end
                
                if closeflag == 1               % check if a close pair has been found
                    removes(1,m) = j;         % index the pair in a n,2 matrix
                    removes(1,m+1) = k;
                    m = m+2;                    % prepare for next pair if found
                end
            end
            
            if xpeak(j)-ceil(sigma2)-1 <= 0 || xpeak(j)+ceil(sigma2)+1 > xw || ypeak(j)-ceil(sigma2)-1 <= 0 || ypeak(j)+ceil(sigma2)+1 > yw
                removes(1,m) = j;
                m=m+1;
            end
        end                                     % When complete all overlapping molecules should be found
        
        if removes(1,1) ~=0                     % triggered if overlapping molecules were found
            ypeak(removes) = [];              % remove row index the list
            xpeak(removes) = [];              % remove column index from list
        end                                     % at this point all overlapping molecules should have been removed
        
        %% Potential localization analysis
        if numel(xpeak)  >0
            for pm = 1:numel(xpeak)
                high_pix_val = i2(ypeak(pm),xpeak(pm));   % grab high pixel value
                i3(:,:) = i2(ypeak(pm)-ceil(sigma2)-1:ypeak(pm)+ceil(sigma2)+1,xpeak(pm)-ceil(sigma2)-1:xpeak(pm)+ceil(sigma2)+1);  % grabs the region around the psf
                %             imagesc(i2);
                %             colormap gray
                %             colorbar
                %             axis image
                %             drawnow
                
                [highrow, highcol] = find(i3(:,:) == max(max(i3(:,:))),1);                  % finds the brightest pixel
                highpix = max(max(i3(:,:)));                                           % records valuei2
                
                
                clear xguess yguess peakguess r0_allguess offguess lowx lowy hix hiy  %clears variables that will be used
                xguess = xgrid(highrow,highcol);                                % initial value for the x position
                yguess = ygrid(highrow,highcol);                                % initial value for the y position
                
                
                % This section is dedicated to guessing the value of N
                lowx = round(xguess-2*ceil(sigma2));    % this region selects
                hix = round(xguess+2*ceil(sigma2));     % bounds on the psf
                lowy = round(yguess-2*ceil(sigma2));    % to sum over to guess
                hiy = round(xguess+2*ceil(sigma2));     % the value of N
                
                % if any values are out of bounds, correct them
                if round(lowx) <=0
                    lowx = 1;
                end
                if round(hix) >= max(max(xgrid))
                    hix = max(max(xgrid));
                end
                if round(lowy) <=0
                    lowy = 1;
                end
                if round(hiy) >= max(max(ygrid))
                    hiy = max(max(ygrid));
                end
                Nguess = sum(sum(i3(lowy:hiy,lowx:hix)));   % initial guess for total number of photons
                
                off_guess = bkgn;
                if bkgn == 0
                    off_guess = 1;
                end
                mol_num = mol_num+1;
                % create an array of input vectors
                beta0(1,mol_num) = xguess;
                beta0(2,mol_num) = yguess;
                beta0(3,mol_num) = Nguess;
                beta0(4,mol_num) = sigma2/2;
                beta0(5,mol_num) = sigma2/2;
                beta0(6,mol_num) = off_guess;
                xpeaks(mol_num,1) = xpeak(pm);
                ypeaks(mol_num,1) = ypeak(pm);
                framenum_temp(mol_num,1) = framenum;
                i4(:,:,mol_num) = i3(:,:);
                clear i3;
                
                if mol_num == 1000
                    
                    [xf_temp, yf_temp, N_temp, off_temp, sigx_temp, sigy_temp,  xf_crlb_temp, yf_crlb_temp, N_crlb_temp, off_crlb_temp, sigx_crlb_temp, sigy_crlb_temp, llv_temp, framenum_all_temp] = chain_loc_widths_crlb(i4, beta0,framenum_temp, xpeaks, ypeaks, xgrid, ygrid, threadcount);
                    
                    clear framenum_temp beta0
                    if  runs ==1
                        xf_all = xf_temp(framenum_all_temp >0);
                        yf_all = yf_temp(framenum_all_temp >0);
                        N = N_temp(framenum_all_temp >0);
                        off_all = off_temp(framenum_all_temp >0);
                        sigx_all = sigx_temp(framenum_all_temp >0);
                        sigy_all = sigy_temp(framenum_all_temp >0);
                        xf_crlb = xf_crlb_temp(framenum_all_temp >0);
                        yf_crlb = yf_crlb_temp(framenum_all_temp >0);
                        N_crlb = N_crlb_temp(framenum_all_temp >0);
                        off_crlb = off_crlb_temp(framenum_all_temp >0);
                        sigx_crlb = sigx_crlb_temp(framenum_all_temp >0);
                        sigy_crlb = sigy_crlb_temp(framenum_all_temp >0);
                        framenum_all = framenum_all_temp(framenum_all_temp >0);
                        llv = llv_temp(framenum_all_temp>0);
                        xc = xpeaks(framenum_all_temp>0);
                        yc = ypeaks(framenum_all_temp>0);
                        i5 = i4(:,:,framenum_all_temp>0);
                        i6 = i4(:,:,framenum_all_temp<0);
                        %                         framenum
                        disp([num2str(100*framenum/(1+n_end-n_start)),'% complete']);
                        runs = runs+1;
                        mol_num = 0;
                        clear i4 xpeaks ypeaks;
                        
                    else
                        xf_all = vertcat(xf_all, xf_temp(framenum_all_temp >0));
                        yf_all = vertcat(yf_all,yf_temp(framenum_all_temp >0));
                        N = vertcat(N,N_temp(framenum_all_temp >0));
                        off_all = vertcat(off_all,off_temp(framenum_all_temp >0));
                        sigx_all = [sigx_all;sigx_temp(framenum_all_temp>0)];
                        sigy_all = [sigy_all;sigy_temp(framenum_all_temp>0)];
                        xf_crlb = vertcat(xf_crlb,xf_crlb_temp(framenum_all_temp >0));
                        yf_crlb = vertcat(yf_crlb,yf_crlb_temp(framenum_all_temp >0));
                        N_crlb = vertcat(N_crlb,N_crlb_temp(framenum_all_temp >0));
                        off_crlb = vertcat(off_crlb,off_crlb_temp(framenum_all_temp >0));
                        sigx_crlb = [sigx_crlb; sigx_crlb_temp(framenum_all_temp >0)];
                        sigy_crlb = [sigy_crlb; sigy_crlb_temp(framenum_all_temp >0)];
                        framenum_all = vertcat(framenum_all, framenum_all_temp(framenum_all_temp >0));
                        llv = vertcat(llv, llv_temp(framenum_all_temp > 0));
                        
                        xc = [xc; xpeaks(framenum_all_temp>0)];
                        yc = [yc; ypeaks(framenum_all_temp>0)];
                        %                         framenum
                        disp([num2str(100*framenum/(1+n_end-n_start)),'% complete']);
                        runs = runs+1;
                        i5 = cat(3,i5,i4(:,:,framenum_all_temp>0));
                        i6 = cat(3,i6,i4(:,:,framenum_all_temp<0));
                        clear i4;
                        mol_num = 0;
                        clear i4 xpeaks ypeaks
                    end
                end
            end
        end
    end
end
if mol_num > 0
    
    [xf_temp, yf_temp, N_temp, off_temp, sigx_temp, sigy_temp,  xf_crlb_temp, yf_crlb_temp, N_crlb_temp, off_crlb_temp, sigx_crlb_temp, sigy_crlb_temp, llv_temp, framenum_all_temp] = chain_loc_widths_crlb(i4, beta0,framenum_temp, xpeaks, ypeaks, xgrid, ygrid, threadcount);
    
    if runs ==1     % Only save values who have a positive framenum all value
        xf_all = xf_temp(framenum_all_temp >0);
        yf_all = yf_temp(framenum_all_temp >0);
        N = N_temp(framenum_all_temp >0);
        off_all = off_temp(framenum_all_temp >0);
        sigx_all = sigx_temp(framenum_all_temp >0);
        sigy_all = sigy_temp(framenum_all_temp >0);
        xf_crlb = xf_crlb_temp(framenum_all_temp >0);
        yf_crlb = yf_crlb_temp(framenum_all_temp >0);
        N_crlb = N_crlb_temp(framenum_all_temp >0);
        off_crlb = off_crlb_temp(framenum_all_temp >0);
        sigx_crlb = sigx_crlb_temp(framenum_all_temp >0);
        sigy_crlb = sigy_crlb_temp(framenum_all_temp >0);
        framenum_all = framenum_all_temp(framenum_all_temp >0);
        llv = llv_temp(framenum_all_temp>0);
        xc = xpeaks(framenum_all_temp>0);
        yc = ypeaks(framenum_all_temp>0);
        %         framenum
        disp([num2str(100*framenum/(1+n_end-n_start)),'% complete']);
        runs = runs+1;
        mol_num = 0;
        i5 = i4(:,:,framenum_all_temp>0);
        i6 = i4(:,:,framenum_all_temp<0);
        clear i4 xpeaks ypeaks
    else
        xf_all = vertcat(xf_all, xf_temp(framenum_all_temp >0));
        yf_all = vertcat(yf_all,yf_temp(framenum_all_temp >0));
        N = vertcat(N,N_temp(framenum_all_temp >0));
        off_all = vertcat(off_all,off_temp(framenum_all_temp >0));                        sigx_all = [sigx_all;sigx_temp(framenum_all_temp>0)];
        sigy_all = [sigy_all;sigy_temp(framenum_all_temp>0)];
        xf_crlb = vertcat(xf_crlb,xf_crlb_temp(framenum_all_temp >0));
        yf_crlb = vertcat(yf_crlb,yf_crlb_temp(framenum_all_temp >0));
        N_crlb = vertcat(N_crlb,N_crlb_temp(framenum_all_temp >0));
        off_crlb = vertcat(off_crlb,off_crlb_temp(framenum_all_temp >0));
        sigx_crlb = [sigx_crlb; sigx_crlb_temp(framenum_all_temp >0)];
        sigy_crlb = [sigy_crlb; sigy_crlb_temp(framenum_all_temp >0)];
        framenum_all = vertcat(framenum_all, framenum_all_temp(framenum_all_temp >0));
        llv = vertcat(llv,llv_temp(framenum_all_temp>0));
        xc = [xc; xpeaks(framenum_all_temp>0)];
        yc = [yc; ypeaks(framenum_all_temp>0)];
        %         framenum
        disp([num2str(100*framenum/(1+n_end-n_start)),'% complete']);
        runs = runs+1;
        i5 = cat(3,i5,i4(:,:,framenum_all_temp>0));
        i6 = cat(3,i6,i4(:,:,framenum_all_temp<0));
        mol_num = 0;
        clear i4 xpeaks ypeaks
    end
end
% subplot(4,4,3); histogram(N); title('Histogram of N')
% subplot(4,4,4); histogram(N_crlb); title('Histogram of N_crlb');
% subplot(4,4,7); histogram(off_all); title('Histogram of Offset');
% subplot(4,4,8); histogram(off_crlb); title('Histogram of Off_crlb');
% subplot(4,4,11); histogram(N_crlb.^0.5./N); title('Fractional Uncertainty in N')
% subplot(4,4,12); histogram(off_crlb.^0.5./off_all); title('Fractional Uncertainty in Offset')
% subplot(4,4,15); histogram(xf_crlb); title('Histogram of Xf crlb')
% subplot(4,4,16); histogram(yf_crlb); title('Histogram of Yf crlb')
% lp2=((sigma2*q/2).^2+(q^2)/12)*1./N+8*pi*((sigma2*q/2).^4)*(bkgn^2)/(q^2)*1./(N.*N);
% lp = 1.3.*sqrt(lp2);
% close(hwait);
clear gpu beta0 button check1 closeflag current_dir distance frame_mol_1 high_pix_val highcol highpix highrow hix hiw i1 i2
clear imagfo j k lowx lowy m max_ims mem_size_per_frame NGuess num_high_pix pm pow removes threadcount xw yw yguess xguess
clear hwait i i4 i_ball i_gauss i_hood ims iprod kd kw llv_temp mem_size mol_num N_temp Nguess off_crlb_temp
clear ypeak ypeaks Ygs ygrid yf_temp yf_crlb_temp xpeaks xpeak Xgs xgrid xf_temp xf_crlb_temp wgrid se rows rk rball off_temp
clear chunkim chunks cols framenum FWHM gpud hiy off_guess runs N_crlb_temp framenum_all_temp framenum_temp
clear a1 amem chunk_size fpath framelist mol_locs max_neur look neur_wait
clear num_chunk num_mol sigma sigma2 sigma2_um tlink

savename = [fname(1:end-4),'_frames_',num2str(n_start),'_',num2str(n_end),'_t_neuro_with_psf.mat'];
save([an_dir,savename]);
end

