%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Learning
%
% This script will generate a random set of neural thetas, analyze an
% image, check the quality of the data, use that to relearn neural thetas,
% and repeat the analysis. This will repeat until a currently undetermined
% point.
%
%
%
% AJN 3-12-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% User Controlled Variables
sp = 50;
lambda = 0.3;
trainper = 0.9;
tossper = 0.05;
epsint = 0.12;
wvlngth = 575;                                   %  peak wavelength of fluorophore emission
NA = 1.4;
q = 0.128;             % Pixel Size
n_start = 1;        % Starting frame
n_end = 0;             % End 0 for end frame set to 0 for final frame
n_bkgn = 100;          % Frame to use for background measurement
pix2pho = 12.78;          % Pixel to photon ratio (for sCMOS set to 1 after pre-processing)
n_start = 1;
n_end = 0;
bkgn = 1;

%% Tolerance Variables
N_tol_min = 10;         % minimum number of photons
N_tol_max = 5000;       % maximum number of photons
off_min = -10;            % minimum number of offset photons
off_max = 12;           % maximum number of offset photons
N_crlb_min = 0;         % minimum variance in number of photons
N_crlb_max = 15000;       % maximum variance in number of photons
sigma_min = 1;        % minimum sigma value
sigma_max = 10;         % maximum sigma value
xf_crlb_min = 0;        % minimum error in number of x-position
xf_crlb_max = 0.91;      % maximum error in number of x-position
yf_crlb_min = 0;        % minimum error in number of y-position
yf_crlb_max = 0.91;      % maximum error in number of y-position
off_crlb_min = 0;     % minimum error in number of offset photons
off_crlb_max = 2;     % maximum error in number of offset photons
sig_crlb_min = 0;
sig_crlb_max = 2;
fr_unc_N = 0.5;           % fractional uncertainty in N
fr_unc_off = 0.5;       % Fractional uncertainty in offset
fr_unc_sig = 0.5;       % Fractional uncertainty in width
x=[];
y=[];

% Code Begins
% build random thetas
theta1 = rand(30, 50)*2*epsint - epsint;
theta2 = rand(1, 31)*2*epsint - epsint;



% load('C:\Users\AJN Lab\Desktop\5-30-17 munc13-halo\Subfolder\Neural_thetas_it_19.mat');
% theta1 = Theta1;
% theta2 = Theta2;
%% In principle we would choose a random file but for now
[fname, fpath] = uigetfile('*.tif');
cd(fpath);
finfo = dir('*.tif');
% save these thetas
% save('neural_thetas_it_0.mat','theta1','theta2');
data_dir = fpath;
an_dir = fpath;
it = 1;
h1 = figure('Units','Normalized','Outerposition',[0 0 1 1]);
x1 = [];
y1 = [];
county = 1;

% amem = 5.9367*10^9; % This is for the Geforce 1060
amem = 956174336; % this is for the geforce 550 Ti
while true
    % analyze file and return fitting variables
    if it < 5
        n_end = 10;
        sp = 2;
    elseif it < 15 && it >=5
        n_end = 50;
        sp = 10;
    elseif it < 25 && it >=15
        n_end = 500;
        sp = 50;
    else
        n_end = 1000;
        sp = 50;
    end
    
    %     figure('Units','Normalized','OuterPosition',[0 0 1 1]);
    [savename] = func_Quhzx_02_7_for_learning(data_dir,finfo(randi(numel(finfo))).name,an_dir,bkgn,theta1, theta2,q,pix2pho,n_start,n_end, wvlngth, NA, amem, sp);
    %     saveas(gcf,['figure_it_', num2str(it),'.tif']);
    %     p1 = pca;
    %     M(county) = getframe(h1);
    %     county = county +1;
    % Tolerancing file
    try
        
        
        %         if it < 10
        %         [x1, y1] = app_gpu_tol_all_color_learning_easy(data_dir, savename, sigma_min, sigma_max, sig_crlb_min, sig_crlb_max, fr_unc_sig, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off);
        %         elseif it >= 10
        [x1, y1] = app_gpu_tol_all_color_learning(it + 20,data_dir, savename, sigma_min, sigma_max, sig_crlb_min, sig_crlb_max, fr_unc_sig, N_crlb_min, N_crlb_max, xf_crlb_min, xf_crlb_max, yf_crlb_min, yf_crlb_max, off_crlb_min, off_crlb_max, N_tol_min,N_tol_max, off_min, off_max, fr_unc_N, fr_unc_off);
        %         end
        %% add tole
        x = [x;x1];
        y = [y;y1];
    catch lsterr
    end
    sum(y)
    %     if sum(y) < 40 && it == 1 % if there are no good fits found, put in images of gaussians
    %         [x,y] = func_build_gauss(x,y, N_tol_min, N_tol_max, sigma_min, sigma_max);
    %     end
    if sum(y) > 40 || it > 1
        
        %% At this point x contains a group of all images and y is a key of whether the image passed the tolerances or not
        %         disp(['Iteration number: ', num2str(it-1), ' and ', num2str(100*sum(y1)/numel(y1)),' % good regions found and ', num2str(sum(y1)), '# of regions found']);
        progress(it,1) = it -1;
        progress(it,2) = sum(y1)/numel(y1);
        progress(it,3) = sum(y1);
        
        
        % now the training set should contain positive and negative examples
        if sum(y)/numel(y) < 1.99 % if neural net picked out less than 70% good regions train new set
            [theta1, theta2] = func_neural_teach( theta1, theta2, x, y, lambda, it, trainper);
            it = it +1;
            
            
            tindex = find(y == 1);
            tosspoint = round(tossper*numel(tindex));
            index = randperm(numel(tindex));
            if it -1 >10
                x(index(1:tosspoint),:) = [];
                y(index(1:tosspoint)) = [];
            end
        else
            break
        end
        %         subplot(1,2,2);
        %         plot(progress(:,1),progress(:,2)*100,'.b','MarkerSize',30);
        %         xlabel('Iteration Number')
        %         ylabel('Percentage of Good Regions')
    else
        theta1 = rand(30, 82)*2*epsint - epsint;
        theta2 = rand(1, 31)*2*epsint - epsint;
%         [x,y] = func_build_gauss(x,y, N_tol_min, N_tol_max, sigma_min, sigma_max);
    end
    
    %     M(it-1) = getframe(gcf);
end