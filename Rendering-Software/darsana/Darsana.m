% Get a good preclean in
close all; clearvars; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Darsana
% Darsana is an auspicious sight of a holy person which bestows merit on
% the person who is seen.
% This is the rendering suite for super res images
% AJN 1/26/18 Ryan Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Variables
% Change these settings to get the image that you want
% Default values are decided by the algorithm and are based off the fitting
% variables for your data set
%%%%% It is assumed that final image touch ups will be done in
%%%%% photoshop
gain = 1;      % gain of final image, effectively changes brightness
smooth = 0.25;   % relative smoothing radius 1 is default value, smooth < 1 = sharpen, smooth > 1 = blur
dz = 50;       % Desired z resolution in nm
grid = 0;       % Grid size, set to 0 for default value
image = 5;      % Image setting variable 1: Grayscale 2:User Color 3:3D 4:User Variable 5: NN render
thresh = 100;   % threshold for NN representation if image is set to 5
vary = 'N';
%% Color Settings
rgb = [255 255 255];  % Color code

% Scale Bar Settings
um_height = 1/(3*1.61803398875); % Height of scale bar in final pixels using the golden ratio
bx = 250;   % xposition of scale bar corner
by = 250;   % y position of scale bar corner, the program will go from [bx, by] -> [bx + um, by + um_height]


% [file, fpath] = uigetfile('*dast*'); % look for a da_storm localized file
files = dir('*dast*');
for mk = 1:numel(files)
mkdir('render');
switch image
    case 1
        [im1, umb] = make_pics(files(mk).name, smooth, grid); % grayscale image function
    case 2
        [im1, umb] = make_colorpics(files(mk).name, rgb/255, smooth, grid);
    case 3
        [im1, umb] = make_3dpics(files(mk).name, dz, smooth, grid); % 3D image function
    case 4
        [im1, umb] = make_varpics(files(mk).name, vary, dz, smooth, grid); % 3D image function
    case 5
        [im1, umb] = make_nnpics(files(mk).name, thresh, dz, smooth, grid); % 3D image function
    otherwise
        im1 = 1;
end
im1 = gain*im1;
%% Create the scale bar
um = round(umb);
mu = round(um_height*um);

% Draw the scale bar
for i = by:by+mu % spanning y space
    for j = bx:bx+um % spanning x space
        im1(i,j,1) = 1;
        if image ~= 1
            im1(i,j,2) = 1;  % set each channel to 1 causing white
            im1(i,j,3) = 1;
        end
    end
end
% while true

if image ~= 1
    subplot(3,4,[1:3;5:7;9:11]);
    imagesc(im1*gain);
    axis image
    % Creating RGB histograms
    imr = im1*0;
    img = imr;
    imb = img;
    imr(:,:,1) = im1(:,:,1)*gain; % red channel
    img(:,:,2) = im1(:,:,2)*gain; % green channel
    imb(:,:,3) = im1(:,:,3)*gain; % blue channel
    subplot(3,4,4);
    imagesc(imr);
    title('Red')
    axis image
    subplot(3,4,8);
    imagesc(img);
    axis image
    title('Green');
    subplot(3,4,12);
    imagesc(imb);
    axis image
    title('Blue');
else
    imagesc(im1*gain);
    colormap('gray');
end


imwrite(uint16(65535*im1),['render\',files(mk).name(1:end-4),'_sm_',num2str(smooth),'_um_',num2str(round(umb)),'.tif']);
end