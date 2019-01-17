%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diff Drift Correction
% A correction algorithm to fix the drift in diffraction limited images
% AJN 7/30/18 Ryan Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; % cleanup line

try
files = dir('*tif'); % grab all tif files (this assumes all files in a folder are apart of the same experiment
for i = 1:numel(files) % load images
    i1(:,:,i) = mean(readtiff(files(i).name),3); % average over all frames in each file
    id = strfind(files(i).name,'_');
    ind(i) = str2num(files(i).name(7:end-4)); % grab index number of files
end
catch lsterr
    [i1, fname, fpath] = readtiff();
    [~,~,o] = size(i1);
    ind = 1:o;
    files = ind;
   
end
%  imagesc(sum(i1,3));
%     [x,y] = ginput(1);
%     x = round(x);
%     y = round(y);
%     wind = -20:20;
%     isub = i1(y+wind, x+ wind,:);
% %     clear i1
%     for i = 1:numel(isub(1,1,:))
%     im1(:,:,i) = interp2(isub(:,:,i),3,'spline');
% %     im1(:,:,i) = imresize(isub(:,:,i),5);
%     end
%     i1 = isub;
%     clear isub
% imagesc(sum(im1,3)) % display 'sum' or uncorrected projection of images
[B,I] = sort(ind); % determine numerical order of files
figure
% for i = 1:numel(files) % create movie of consecutive frames
    i2 = i1(:,:,I(:)); % build array i2 to be a chronological array of images
%     % display stuff and grab frame
%     imagesc(i2(:,:,i))
%     drawnow
%     title('Uncorrected Image');
%     xlabel(['Frame number ',num2str(i)])
%     M(i) = getframe(gcf);
% end
% movie2gif(M,'movie.gif','DelayTime',0.05,'LoopCount',Inf)
[drifts] = get_drift_ims(i2); % fix the drift of the images
% figure
% [m1, n1, o1] = size(im1);
% [m2, n2, o2] = size(isub);
drifts = drifts;
% smooth drifts w/ gaussian filter
drifts(:,1) = gausssmooth(drifts(:,1),4,10);
drifts(:,2) = gausssmooth(drifts(:,2),4,10);
plot(drifts(:,1),drifts(:,2)) % plot drifts in x,y
figure
% save([fname(1:end-4),'_drifts.mat'],'drifts');
% approximate correction with circshift function and display result
for i = 1:numel(files)
    subplot(1,2,1);
    imagesc(i2(:,:,i))
    axis image
    title('Uncorrected')
    i3(:,:,i) = circshift(imresize(i2(:,:,i),10),-round(10*[drifts(i,2),drifts(i,1)])); % we resize the image to reduce the approximation that comes in with rounding
    subplot(1,2,2);
    imagesc(i3(:,:,i));
    title('Corrected');
   
    axis image
    drawnow
     N(i) = getframe(gcf);
    
end