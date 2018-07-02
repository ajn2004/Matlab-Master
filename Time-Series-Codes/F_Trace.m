clearvars; clc; close all;
%% F_trace
% This script will allow the user to select a number of points and create a
% series of images showing the average F over a user selected number of
% boutons as well as a where cutoffs are to measure the noise associated
% with that signal. This should work with any image file type

%% USER VARIABLES
im_type = 'tif';  % image extension currently either 'tif' or 'fits'
stim_fr = 100;        % First frame stimulation occurs on
stims = 1;            % Total number of stimulations
str = 10;            % Stimulation rate in Hz
pixw = 7;           % Pixel width in radius (i.e. value of 7 gives 15x15 square window)
%%%%%%%%%%%%%%%%%%%%%END USER CONTROL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analysis
c = scrub_config();  % reads camera configuration file
[i1, mname, mpath] = read_im(im_type); % loads in an image and corrects
[m,n,o] = size(i1); % get image size in [rows, cols, frames]

tex = c.AccumulateCycleTime;

si1 = std(i1,1,3);  % create standard deviation image for selection
imagesc(si1)    % display standard deviation image
axis image

if exist([mpath,mname(1:end-4),'_ROIlist.mat'])
    roi_load = input('A ROI file has been found, would you like to load it(y/n)?', 's');
    if strcmp(roi_load,'y') || strcmp(roi_load,'Y')
        load([mpath,mname(1:end-4),'_ROIlist.mat']);
    end
end

if ~exist('x')
    while true  % ensure user must input a number
        try
            num = input('How many points?');  % ask user how many points to select
            break
        catch lsterr
        end
    end
    
    for i = 1:num  % Select regions and create boxes around selections
        [x(i),y(i)] = ginput(1);
        draw_boxes([x(i),y(i)],pixw);
    end
    
    wind = -pixw:pixw;  % define a variable for easy sub_image selection
    
    for i = 1:num  % select out each sub region and measurements
        sub1 = i1(round(y(i)) + wind, round(x(i)) + wind,:);  % grab image subregion
        sfl = sum(sum(sub1));    % sum up region fluorescence in each frame
        sfluor = sum(sum(sub1)); % creates a 1x1xo variable of total fluorescence/region/frame
        sfluor = sfluor(:);   % linearize array
        ifluor(:,i) = sfl(:);  % ifluor is for individual fluor
    end
    mfluor = mean(ifluor,2);  % average over all individuals gives mfluor (mean fluorescence)
end
t = (1:o)*tex;

%% DATA PLOTTING
plot(t,ifluor,'color', [0, 0, 0] + 0.75); % make individuals traces gray
hold on
plot(t,mfluor, 'color', 'k');  % make mean trace black
% plot stims
for i = 1:stims
    plot([stim_fr*tex + (i-1)/str,stim_fr*tex + (i-1)/str],[min(mfluor),max(mfluor)],'r'); % plot a red line at stim_frame * s/frame + (stimnumber-1)/stimspersec
end
hold off
xlabel('Time in [s]')
ylabel('F in [A.U.]');
title('Trace of mean F with individuals');
% create new figure of just average fluorescence
figure % make new figure
plot(t,mfluor,'k');  % plot only average trace
hold on
% plot stims
for i = 1:stims
    plot([stim_fr*tex + (i-1)/str,stim_fr*tex + (i-1)/str],[min(mfluor),max(mfluor)],'r'); % plot a red line at stim_frame * s/frame + (stimnumber-1)/stimspersec
end
hold off
xlabel('Time in [s]')
ylabel('F in [A.U.]');
title('Trace of mean F');


save([mpath,mname(1:end-4),'_ROIlist.mat'],'x','y','ifluor','mfluor');

