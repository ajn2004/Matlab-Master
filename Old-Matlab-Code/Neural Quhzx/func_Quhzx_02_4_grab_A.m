function func_Quhzx_02_4_grab_A(fpath,fname,an_dir,bkgn,Theta1, Theta2, q,pix2pho,n_start,n_end, wvlngth, NA, n_bkgn)
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
imagfo = imfinfo([fpath,fname]);                                          % variable to measure image info
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
sigma2_um = 0.55*wvlngth/(1000*NA);         % width of the Rayleigh radius in um
sigma2 = sigma2_um/q;                       % width in pixels
sigma = sigma2/2;                           % standard deviation of the gaussian
[xgrid, ygrid] = meshgrid(-ceil(sigma2)-1:ceil(sigma2)+1,-ceil(sigma2)-1:ceil(sigma2)+1);   % grids used for building fitting gaussians
wgrid = 2*ceil(sigma2)+3;                                                                   % width of the grids for fitting


runs = 1;

threadcount = 1024;
if n_end == 0
    n_end = numel(imagfo);
end


i1 = zeros(imagfo(1).Height,imagfo(1).Width,n_end-n_start+1);
for i = n_start:n_end
    j = i-n_start +1;
    i1(:,:,j) = double(imread([fpath,fname],i,'Info',imagfo))/pix2pho;
end
disp('Done Loading');
% disp('Past Loading');
% tic
%% Background Subtraction
rball=5; %radius of rolling ball

se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball
i_ball = se.getheight();
i_hood = se.getnhood();
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
FWHM = 1;
rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
% xrow = -kw/2:kw/2;
i_ball= i_ball./sum(sum(i_ball));

% gauss_vec = exp(-2*xrow.*xrow/(rk*rk));
i_gauss=exp(-2*kd.*kd./(rk*rk));
i_gauss=i_gauss/sum(sum(i_gauss));

% calculate memory requirements
gpud = gpuDevice;
rows = numel(i1(:,1,1));
cols = numel(i1(1,:,1));
ims = numel(i1(1,1,:));

mem_size = rows*cols*ims*8;

if mem_size < gpud.AvailableMemory / 8             % if memory of the image is small enough, run it through
    [iprod] = image_process(i1,i_gauss, i_ball);
else                                                % If memory of image is too small, break it down into chunks and run it through the gpu
    chunks = ceil(mem_size / (gpud.AvailableMemory/10));
    chunkim = ceil(ims / chunks);
    for i = 1:chunks
        if i ~= chunks
            [iprod_temp] = image_process(i1(:,:,1+(i-1)*chunkim:i*chunkim),i_gauss, i_ball);
            %             [atemp] = image_neural_2(iprod_temp, Theta1.', Theta2.', numel(iprod_temp(1,1,:)));
            if i ==1
                iprod = iprod_temp;
                %                 a1 = atemp;
                clear iprod_temp
            else
                iprod = cat(3,iprod, iprod_temp);
                %                 a1 = cat(3,a1,atemp);
                clear iprod_temp atemp
            end
        else
            [iprod_temp] = image_process(i1(:,:,1+(i-1)*chunkim:end),i_gauss, i_ball);
            %             [atemp] = image_neural_2(iprod_temp, Theta1.', Theta2.', numel(iprod_temp(1,1,:)));
            iprod = cat(3,iprod,iprod_temp);
            %             a1 = cat(3,a1,atemp);
            clear iprod_temp
        end
    end
end
disp('Done Subtracting background');

max_neur = 135*135*200*8; %this number is based off simulations

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
        [atemp] = image_neural_2(iprod(:,:,1+ (i-1)*chunk_size:end), Theta1.', Theta2.', numel(iprod(1,1,1+ (i-1)*chunk_size:end)));
        a1 = cat(3,a1,atemp);
        clear atemp
    end
    waitbar(i/num_chunk,neur_wait,['Neural Net calcualtion is ', num2str(100*i/num_chunk),'% complete']);
end
close(neur_wait);
disp('Neural calc done');
clear i1

clearvars -except a1 an_dir fname iprod

save([an_dir,fname(1:end-4),'_activation.mat'],'a1', 'iprod');
end

