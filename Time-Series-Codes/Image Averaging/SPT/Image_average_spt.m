% Image Averager
%
% A simple script to take in several time-averaged trials and create 1
% final averaged copy
%
% Ajn 3/7/17
%
% Version 2 will automatically combine all files in a folder and
% additionally output a single image of the peak frame - the average frames
startitup;
clearvars

psize = 6;
rball = 5;
stim = 25;
wvlngth = 509;
NA = 1.4;
q = 0.156;
pix2pho = 33.6;




[fname, fpath] = uigetfile('*.tif','Choose a baseline tiff'); % Allow user to select a director with images in it

cd(fpath); % change working directory to chosen directory
files = dir('*.tif'); % get a list of all images in directory

sigma2_um = 0.61*wvlngth/(1000*NA);         % width of the Rayleigh radius in um
sigma2 = sigma2_um/q;                       % width in pixels
sigma = sigma2/2;                           % standard deviation of the gaussian
for j = 1:numel(files) % loop over every image found
    finfo = dir([files(j).name(1:end-4),'*']); % find images with the same base file name, as andor labels _1, _2, ...
%     if numel(finfo) > 1
%         for i = 1:numel(finfo)
            %% Image Averaging
            i1 = readtiff(files(j).name);
            %     i1 = readtiff(fname);
            
            [iprod] = func_background_sub(i1/pix2pho, 3*sigma2, rball);
            ave2 = mean(iprod(:,:,1:stim-1),3);
            [ave3] = ave_subtract(iprod, stim);
%             [theta1, theta2] = func_minutia(ave3, sigma2,files(j).name);
            %     ave2 = mean(i1(:,:,1:stim-1),3);
            [xf_all, yf_all, N_all, sigx_all, sigy_all, off_all, framenum_all, number, xf_crlb, yf_crlb, N_crlb, sigx_crlb, sigy_crlb, off_crlb,llv, points] = func_time_series_spt(ave3, psize, sigma2);
            save([files(j).name(1:end-4),'_results.mat'],'xf_all', 'yf_all', 'N_all', 'sigx_all', 'sigy_all', 'off_all', 'framenum_all', 'number', 'xf_crlb', 'yf_crlb', 'N_crlb', 'sigx_crlb', 'sigy_crlb', 'off_crlb','llv', 'points');
        end
%     end
% end