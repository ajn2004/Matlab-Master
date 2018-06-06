%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region Measurer
% Allows user to crop a region then create a bit map based off a density
% plot to allow matlab to use the region_props function
%
%
% 9/17/15 AJN LW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% numovert = 15;  %number of verticies used to select a cell
% grid_size = 30; % in nm
q =1 ;%pixel size in nm

[fname, fpath] = uigetfile('*.tif','Select a rendered file');  % grab file to be analyzed
cd(fpath);
i2 = imread([fpath,fname]); %load desired file

i1 = rgb2gray(i2);
%% Show the data
% f1 = figure;
% screen_size = get(0, 'ScreenSize');
% set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% h1 = imagesc(i1, [min(i1(:)) max(i1(:))]);
% axis equal
% colormap('gray');
% title('Zoom to preference, press enter when ready');
% p=gca;  % get current axes
% 
% 
% hold on
% % maskx = zeros(size(
% w = 1; % counting variable
% result = input('Press enter when ready');
% %% creat a polygon that encloses the roi
% 
% %% creat a polygon that encloses the roi
% title('Select the verticies of the region of interest for polygon, press enter when complete');
% % continue to select points until a click is placed within 100 nm of
% % original point
% while true
%     clearvars points;
%     % selects point clicked in plot
%     bc = waitforbuttonpress;
%     if bc == 1
%         break
%         
%     else
%         points = get(p,'currentpoint');
%         
%         
%         % assigns selected points to array
%         poly_vert1(w,1) = points(1,1);
%         poly_vert1(w,2) = points(1,2);
%         
%         % draws polygon on plot
%         plot(poly_vert1(:,1),poly_vert1(:,2),'r');
%         w = w+1;
%     end
% end
% 
% saveas(f1,strcat(['Selection Mask for ' base_name(1:end-4)]), 'tif');
% 
% 
% poly_index = inpolygon(xf_all*q,yf_all*q,poly_vert1(:,1),poly_vert1(:,2)); % creat an index of molecules inside region defined
% xf_in = xf_all(poly_index)*q;
% yf_in = yf_all(poly_index)*q;
% pol_area = polyarea(poly_vert1(:,1),poly_vert1(:,2));  % measure area of selected polygon
% avg_den = numel(xf_in)/pol_area; % in molecules / square micron
% disp(['Average density is ', num2str(avg_den),'molecules/um^2']); %display function
% avg_den_nm = avg_den/1000^2;
% avg_mol = avg_den_nm * grid_size*grid_size;
% 
% %% Allow user to select grid size based on molecular density
% try
%     while true
%         s = input(['Your grid size is ', num2str(grid_size), 'nm which gives ', num2str(avg_mol),' molecules per grid on average, use this value? y/n: '], 's');
%         if s == 'y' || s == 'Y'
%             break
%         else
%             grid_size = input('Enter new grid size: ');
%             avg_mol = avg_den_nm * grid_size*grid_size;
%         end
%     end
% catch lasterr
% end

%% Divide molecules into grids
% xmax = grid_size*ceil(max(xf_in)*1000/grid_size);
% ymax = grid_size*ceil(max(yf_in)*1000/grid_size);
%build a mesh grid for finding molecule location
% [Xgrid, Ygrid] = meshgrid(0:grid_size:xmax, 0:grid_size:ymax);

% dens = zeros(ymax/grid_size,xmax/grid_size);  % allocate space

clear i  % free counting variable

%% Density filtering
% this is a two step process where the grid size is doubled from chosen
% value and a density filter is performed. The grid size is then haved and
% the density filter is performed again on the remaining molecules
% build density matrix over all molecules in previously drawn polygon
% for i = 1:numel(xf_in)
%     x_ind = find(Xgrid(1,:) > xf_in(i)*1000, 1, 'first') - 1;
%     y_ind = find(Ygrid(:,1) > yf_in(i)*1000, 1, 'first') - 1;
%     dens(y_ind,x_ind) = dens(y_ind,x_ind) + 1;
% end

%% Density filtering by threshold
threshold = 5;
[row, col] = find(i1 < threshold);   %find molecules below threshold
filt_dens = i1;

for j = 1:numel(row)  % I don't know why, but without this loop program freezes
    filt_dens(row(j),col(j)) = 0; % loop over indecies and set to 0
end
imagesc(filt_dens); % show image
title(['Density image filtered at ', num2str(threshold)]);
axis equal

%% Allow user to select density threshold
try
    while true %repeat until user chooses to move on
        s = input(['Your threshold is ', num2str(threshold), ', use this value? y/n: '], 's');
        if s == 'y' || s == 'Y'
            break
        else %if user selects to chose new threshold rerun at new threshold
            threshold = input('Enter new threshold: ');
            [row, col] = find(i1 < threshold);
            filt_dens = i1;
            for j = 1:numel(row)
                filt_dens(row(j),col(j)) = 0;
            end
            dens_filt = medfilt2(filt_dens,[5,5]);
            [row, col] = find(dens_filt < threshold);
            for j = 1:numel(row)
                dens_filt(row(j),col(j)) = 0;
            end
            imagesc(dens_filt);
            axis equal
            title(['Density image filtered at ', num2str(threshold)]);
        end
    end
catch lasterr
end


%% Build Binary image for processing
[rowf, colf] = find(dens_filt > 0); %find values above threshold

bits_Fin = filt_dens.*0; % allocate array

for k = 1:numel(rowf)  %build bitmap
    bits_Fin(rowf(k),colf(k)) = 1;
end

BW = logical(bits_Fin); %convert filtered image to logical for region props

%% Region props and results
s = regionprops(BW, 'all'); % run region props
Centroids = cat(1,s.Centroid);      %
Areas = cat(1, s.Area);             % these lines organize structure data
Diams = cat(1, s.EquivDiameter);    %

area_thresh = 3; %initialize area threshold

%display results of threshold and begin user selection process
imshow(BW);
axis equal
title('Circled regions are being selected for with current threshold');
area_ind = find(Areas > area_thresh);
some_cents = Centroids(area_ind,:);
hold on
plot(some_cents(:,1), some_cents(:,2), 'b.');
radii = Diams(area_ind)./2;
circles(some_cents(:,1),some_cents(:,2), radii , 'r');
hold off
%% Allow user to select area threshold to select against high pixels
try
    while true
        % the following 2 lines are User information and ask for input
        disp(['Your threshold area of ', num2str(area_thresh),' pixels^2 corresponds to ', num2str(area_thresh*q^2/1000000),'um^2']);
        as = input(['Your threshold area is ', num2str(area_thresh), ' pixels^2, use this value? y/n: '], 's');
        if as == 'y' || as == 'Y' %check user response and verify selection
            conf = input(['Use ', num2str(area_thresh), ' pixels^2 as a threshold for area? (y/n): '], 's');
            if as == 'y' || as == 'Y'
                break
            end
        else
            area_thresh = input('Enter new threshold: '); %ask user to choose new threshold
            % update and display results with circles around relevant structures
            imshow(BW);
            axis equal
            area_ind = find(Areas > area_thresh);
            some_cents = Centroids(area_ind,:);
            radii = Diams(area_ind)./2;
            hold on
            plot(some_cents(:,1), some_cents(:,2), 'b.');
            circles(some_cents(:,1),some_cents(:,2), radii , 'r');
            hold off
        end
    end
catch lasterr
end

%% cleanup and save
fields = fieldnames(s); % get names of fields in s struct
for i = 1:numel(fields); % this loop creates variables for all single vector fields in s
    v = genvarname(strcat(fields(i),'_thresh')); %generate variable name based on relevant field name
    try
        eval(char(strcat(v, ' = cat(1, s(area_ind).(fields{i}));'))); %forces matlab to evaluate this string
    catch lasterr
    end
end


save_name = [base_name(1:end-4),'_regions.mat']; %define file name based off previous name

% CLEAN UP VARIABLE LIST
clear Areas Centroids Diams Xgrid Ygrid a0_xcorr an_dir an_dir2 array_size_x array_size_y
clear axes2 axes3 axes4 axes5 base_name base_name2_an base_name2_in base_name2_out
clear checkbox1 col colf corr_unit_mol correct_drift corres_time data_dir
clear data_dir2 data_dir2_an debug_button edit1 edit2 edit3 edit4 edit5 edit6 edit7 edit 8
clear edit9 eval_px_frm eval_px_t eval_py_t eval_py_t eventdata f1 fields figure1 files
clear fname fname_temp fpath fpath_temp fr_un frames_last_file grab grab_sum_all h1 hObject
clear handles clear i imfiles input1 input2 intens j jk k kj lasterr mat_file max_unc mse
clear new_name nfiles num_corr_frames numovert off_xcorr out_file output p
clear polyfitdone proj pushbutton1 pushbutton12 pushbutton13 puhbutton14 pushbutton4 pushbutton9
clear r0_xcorr radii radiobutton3 radiobutton4 radiobutton5 radiobutton6 row rof
clear save_mat_file screen_size shift size_grid size_grid_um slider1 slider2 slider3 some_cents t text10
clear text10 text11 text12 text13 text14 text16 text18 text19 text2 text20 text21 text22 text3 text4 text5
clear text6 text7 text8 threshold tot_name uipanel1 uipanel3 uipanel4  uipanel5 uipanel6 uipanel7 uipanel8 uipanel9
clear v w x0_box x1_box x_ind x_max xcorr_array xcorr_j xcorr_i xdrift_pix xfit_poly_order xind xmax
clear xpix xpix2 y0_box y1_box y_ind y_max ydrift_pix yfit_poly_order ymax ypix ypix2 edit8 eval_py_frm
clear conf as pushbutton14
clc
save(save_name);
close all
