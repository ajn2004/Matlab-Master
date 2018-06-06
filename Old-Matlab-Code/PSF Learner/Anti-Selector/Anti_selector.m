%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anti-Selector
%
% This script will allow the user to select a folder containing activation
% images. The script will load the file and present regions identified by
% the activation image. The user will then be prompted to click regions
% that are likely not molecules. These regions will be recorded and stored
% as a separate training set to be used in subsequent training of the
% neural net.
%
% AJN 3-11-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Up
close all
% clear all
% clc

% Selected Variables
num_frames = 200;
width = 3; 
% Folder Identifications
% 
% [fname, fpath] = uigetfile('*.mat','Select an activation image');
% % [fname, fpathi] = uigetfile('*.tif','Select corresponding Tiff Images');
% 
% % Directory Building
% finfo = dir([fpath,'*.mat']);
% x = [];
% y = [];
for i = 3 : numel(finfo)

   figure('units','normalized','outerposition',[0 0 1 1]);
   load([fpath,finfo(i).name]);
   frames = randperm(numel(iprod(1,1,:)));

   for k = 1: num_frames
       [mol_locs, num_mol] = neural_gpu(a1(:,:,frames(k)));  % identify regions that would be considered

       imagesc(iprod(:,:,frames(k)));
       colormap('gray')
       axis image
       boxes(mol_locs,width+.5,'g'); % draw green boxes around identifications
       pl =  gca; 
       w = 1;
       
       %% User Selects Points to change
       while true  % Allow user to select points
           clearvars points;
%            selects point clicked in plot
           title('Select a new center, press enter to quit')
           bc = waitforbuttonpress;
           if bc == 1 % if user hits enter break the loop
               w = w-1;
               break
           else
               points = get(pl,'currentpoint');
%                assigns selected points to array
               cents(w,1) = points(1,1); % save value of center in um
               cents(w,2) = points(1,2); % save value of center in um
               % Find closest molecule location consistent with the point
               rs = ((mol_locs(:,1) - cents(w,1)).^2 + (mol_locs(:,2)-cents(w,2)).^2).^0.5;
               ind(w) = find(rs == min(rs));
               boxes(mol_locs(ind(w),:),width+.5,'r') % draw a red box around that point
               w = w+1;
           end
       end
       %% Program stores new points and saves negative image
       if w > 0 
           for l = 1:w
               imt = iprod(mol_locs(ind(w),2) - width : mol_locs(ind(w),2) + width,mol_locs(ind(w),1) - width : mol_locs(ind(w),1) + width,frames(k));
               x = [x;imt(:).'];
               y = [y;0];
           end
       end
   end
end
   save('Negative_set.mat','x','y');