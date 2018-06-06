%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A program to measure the pixel to photon ratio of the sCMOS camera
% 7/6/13 AJN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE TO USER: This program assumes that you have taken several images of
% a homegenous light source at different intensities and have labeled them
% sequentially such as example_1.tif. 
% This program will automatically load all files in sequence and analyze them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the workspace
clear all
close all
clc

numochan = 1;
%% Declare some variables
file_num = 6;  % number of files in your series
frame = 100;    % number of frames in each picture
num_of_vert = 5;
%% Load primary file and creat series setup
[fname fpath] = uigetfile('*.tif', ' Choose the file that starts the series');
cd(fpath);
an_dir = 'L:\Data\2-8-14 2color smcos\tIFFS\pix2pho\';
sname = fname(1:numel(fname)-5);        %creates a name for all series
im1 = imread((strcat(fpath,fname)),1);
clims = [min(im1(:)) max(im1(:))];
%% plot the data
f1 = figure;
screen_size = get(0, 'ScreenSize');
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
h=imshow(im1,clims) ;
axis equal
title('Zoom to preference, press enter when ready');
p=gca;  % get current axes

hold on
% maskx = zeros(size(

result = input('Press enter when ready');
%% creat a polygon that encloses the roi
colo = ['r' , 'g'];
for lk = 1:numochan
    title(['Select the verticies of the region of interest for polygon_',num2str(lk)]);
    for w = 1:num_of_vert
        clearvars points;
        % selects point clicked in plot
        waitforbuttonpress;
        points = get(p,'currentpoint');
        
        if lk ==1
        % assigns selected points to array
        poly_vert1(w,1) = round(points(1,1));       
        poly_vert1(w,2) = round(points(1,2));
        
        % draws polygon on plot in red
        plot(poly_vert1(:,1,lk),poly_vert1(:,2,lk),colo(lk));
        elseif lk == 2
            poly_vert2(w,1) = round(points(1,1)); 
            poly_vert2(w,2) = round(points(1,2));
            plot(poly_vert2(:,1),poly_vert2(:,2),colo(lk));
        end
    end
end
saveas(f1,strcat(an_dir,'Selection Mask'), 'tif');
h = waitbar(0,'Loading images 0%');

%% Creation of mean pixel map and average pixel map
for k = 0:file_num-1
    clearvars A
    %% Load all files in series
    for i = 1:frame
        A(:,:,i) = imread(strcat(fpath,sname,num2str(k),'.tif'),i);
        waitbar((k*frame+i)/((file_num)*frame),h,[num2str(100*(k*frame+i)/((file_num)*frame)), '% of data is loaded']);
    end
    %% Creation of average and variance of each pixel
    Aave(:,:,k+1) = mean(double(A),3);       % average
    Avar(:,:,k+1) = var(double(A),0,3);      % variance

end
close(h);
f = waitbar(0,'0% of data analyzed');
c=1;
%% Calculation of pixel to photon ratio
for m = 1:length(A(:,1,1))
    for n = 1:length(A(1,:,1))
        clear p
        p = polyfit(Aave(m,n,:),Avar(m,n,:),1); % uses polyfit to measure slope of mean vs. variance
        pixmap(m,n) = p(1);                     % generates a pixel map of slopes (pix2pho ratio)
        offset(m,n) = p(2);
        waitbar(c/(length(A(:,1,1))*length(A(1,:,1))),f,[num2str(100*c/(length(A(:,1,1))*length(A(1,:,1)))), '% of data is analyzed']);
        c=c+1;
    end
end
avepix = mean(mean(pixmap(find(1-isnan(pixmap)))));
close(f)

%% Setting outside ROI values
[X1, Y1] = meshgrid(1:numel(im1(1,:)),1:numel(im1(:,1)));
[X2, Y2] = meshgrid(1:numel(im1(1,:)),1:numel(im1(:,1)));
INdex1 = inpolygon(X1,Y1,poly_vert1(:,1),poly_vert1(:,2));
INdex2 = inpolygon(X2,Y2,poly_vert2(:,1),poly_vert2(:,2));
INdex = INdex1+INdex2;
pixmap(1-INdex) = 10^5;
save(strcat(fpath,'pix2photo_map'));
