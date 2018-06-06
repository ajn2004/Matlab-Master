%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to estimate the micron per pixel value given a scale image
% 
% This program accepts a tiff image of a scale image. The image should be
% oriented so that the scale bars are lying horizontal, and separated by a
% vertical distance, relative to this screen. (i.e. should look like a
% ladder)
%
% You will be asked questions relating to your optical setup. The code is
% already setup to handle both the Hamamtsu sCMOS and the Andor iXon
% cameras available in the lab.
%
% You are required to put in an expected pixel value based off of a manual
% measurement. This tells the program to look in a window around the
% expected value, and helps reduce errors in the analysis.
%
% One the image is loaded and the questions are answered, you are asked to
% select a region in the image to analyze. This is done by zooming into the
% image allowing you to select the appropriate region for analysis. Once
% zoomed in, press 'enter' and click 5 points (4 points of a rectangle and
% the 5th to enclose the region in a polygon) to indicate which area you
% want analyzed
%
% In my experience choosing a narrow region with many bars works best.
%
% This program can still analyze out of focus images.
%
% AJN 5/4/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Load the data file into the workspace
[fname1, fpath1] = uigetfile('*.tif','Choose a file to analyze');
uncropped = double(imread((strcat(fpath1,fname1))));
addpath(pwd);
cd(fpath1);
clims = [min(uncropped(:)), max(uncropped(:))];
%% plot the data
f1 = figure;
screen_size = get(0, 'ScreenSize');
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
h = imshow(uncropped,clims);
axis equal

scmos = input('Did you use the sCMOS camera for this experiment? y or n ','s');
if strcmp(scmos,'y')||strcmp(scmos,'Y')
    real_pix = 6.4; % real pixel size in microns
else
    real_pix = 16;
end
obj_mag = input('What magnification objective did you use? ');
tel_str = input('What was the approximate magnification of your telescope? ');
expected_pix = input('What is your expected pixel size in nm? ');
totalmag = tel_str*obj_mag;
expect_pix = real_pix/totalmag;
scale_sep = 10; % scale separation size in microns


% get current axes
p=gca;
hold on
% maskx = zeros(size(
%defines the number of vertices one will be selecting
num_of_vert = 5;
result = input('Press enter once zoomed in');
%% creat a polygon that encloses the roi
for w = 1:num_of_vert
    clearvars points;
    % selects point clicked in plot
    disp('Choose an area to measure');
    waitforbuttonpress;
    points = get(p,'currentpoint');
    
    % assigns selected points to array
    poly_vert(w,1) = points(1,1);       
    poly_vert(w,2) = points(1,2);
    
    % draws polygon on plot in red
    plot(poly_vert(:,1),poly_vert(:,2),'r');
    
end
saveas(f1,strcat(fpath1,fname1(1:numel(fname1)-4),'_figure.fig'),'fig');
hold off
% poly_vert = ([6.041099636	11.3361058; 6.364811368	11.54877688 ; 7.731176006	10.07513568; 7.482746072	9.911397765; 6.044863726	11.32481353]);

cropped = uncropped(min(poly_vert(:,2)):max(poly_vert(:,2)),min(poly_vert(:,1)):max(poly_vert(:,1)));

% fst_crop = gpuArray(cropped);
C1_gpu = xcorr2(cropped);
C1 = gather(C1_gpu);
[yind,xind] = find(C1==max(C1(:)));
plot(1:numel(C1(:,1)),C1(:,xind));
title('Autocorrelation along the vertical axis');
ylabel('Correlation a.u.');
xlabel('Displacement in pixels');
center = (numel(C1(:,1))+1)/2;
non_rot = horzcat((1:numel(C1(:,1))).',C1(:,xind));

% Define each point as a vector and split the distribution in half about
% the center
non_rot1 = non_rot(1:center,:);
% non_rot2 = fliplr(non_rot(center:numel(C1(:,1)),:).').';
% non_rot2(:,1) = non_rot1(:,1);

% rotate the vectors by the angle created with the origin and maximum centeral point
theta1 = atan(non_rot1(end,2)/non_rot1(end,1));
rot1 = zeros(numel(non_rot1(:,1)),2);
rot2 = zeros(numel(non_rot1(:,1)),2);
rot1(:,1) = cos(-theta1).*non_rot1(:,1) - sin(-theta1).*non_rot1(:,2);
rot1(:,2) = sin(-theta1).*non_rot1(:,1) + cos(-theta1).*non_rot1(:,2);
figure
plot(rot1(:,1),rot1(:,2))
title('Rotated autocorrelation curve');

% Find the peaks in the distribution
peaks1 = findpeaks(rot1(:,2));
index1 = ismember(rot1(:,2),peaks1);


%% This section is reduntant considering the symmetry of the data
% % do it again for the second half of the distribution
% theta2 = atan(non_rot1(end,2)/non_rot1(end,1));
% 
% rot2(:,1) = cos(-theta2).*non_rot2(:,1) - sin(-theta2).*non_rot2(:,2);
% rot2(:,2) = sin(-theta2).*non_rot2(:,1) + cos(-theta2).*non_rot2(:,2);
% plot(rot2(:,1),rot2(:,2))
% peaks2 = findpeaks(rot2(:,2));
% index2 = ismember(rot2(:,2),peaks2);
% % program's guesses in pixels
guess1 = center-(cos(theta1).*rot1(index1,1)-sin(theta1).*rot1(index1,2));
% guess2 = -center+cos(theta2).*rot2(index2,1)-sin(theta2).*rot2(index2,2);

% all guesses combined in um
allguess = scale_sep./guess1;
% second_guess = scale_sep./guess1(index2);

    [val ind_gue] = min(abs(allguess*1000-expected_pix));

if abs(allguess(ind_gue)*1000-expected_pix)/expected_pix > 0.2;
    ind_gue = find(allguess > 0.8*expect_pix & allguess<1.2*expect_pix );
end
end_guess = allguess(ind_gue);
per_err = abs(end_guess-expect_pix)/expect_pix;
actual_mag = real_pix/end_guess;
disp(['Measured a pixel size of ', num2str(end_guess), 'um which is in ', num2str(per_err*100),'% error of the calculated value']);
disp(['Calculated pixel size: ', num2str(expect_pix)]);
disp(['Expected Magnification: ', num2str(totalmag),'X']);
disp(['Measured Magnifiation: ', num2str(actual_mag),'X']);
%attempts to rotate autocorrelation to find 

% clearvars -except uncropped cropped poly_vert fname1 fpath1 C1
%   
% 
% save(strcat(fname1(1:numel(fname1)-4),'_analyzed.mat'));
% 
% % answer = input('Load a rendered image of total cell? (y/n)','s');
% % if answer == 'y'
% %     [fname2, fpath2] = uigetfile('*.tif','Choose a rendered image');
% %     image1 = imread(strcat(fpath2,fname2));
% %     mask = roipoly(image1,poly_vert(:,1),poly_vert(:,2));
% end
