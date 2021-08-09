function cdata_out = icp_cdata_drift_correct(cdata, mol_chunk, color, max_iter)
cdata_out = cdata;
% Handle optional variables
if ~exist('mol_chunk','var')
    mol_chunk = 100;
end
if ~exist('color','var')
    color = 'red';
end
if ~exists('max_iter','var')
    max_iter = 10000;
end

% Coordinate centering
try % Try orange red combo
    [xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
    xfo = xfo - mean(cdata.red.xf);
    yfo = yfo - mean(cdata.red.yf);
    zfo = zfo - mean(cdata.red.zf);
    zfo = cdata.orange.zf;
    xfr = cdata.red.xf - mean(cdata.red.xf);
    yfr = cdata.red.yf - mean(cdata.red.yf);
    zfr = cdata.red.zf - mean(cdata.red.zf);
catch % If failure, it could be either a lack of red or a lack of orange
    try %Try a red centering
        xfr = cdata.red.xf - mean(cdata.red.xf);
        yfr = cdata.red.yf - mean(cdata.red.yf);
        zfr = cdata.red.zf - mean(cdata.red.zf);
    catch %if that fails, only orange exists
        [xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
        xfo = xfo - mean(xfo);
        yfo = yfo - mean(yfo);
        zfo = zfo - mean(zfo);
    end
end

switch color % color specifies which group of localizations will be used to drift correct
    
    case 'red'
        unsort_locs = [xfr, yfr, zfr, cdata.red.framenumber];
    case 'orange'
        unsort_locs = [xfo, yfo, zfo, cdata.orange.framenumber];
    case 'both'
        unsort_locs = [xfr, yfr, zfr, cdata.red.framenumber];
        unsort_locs = [unsort_locs;[xfo, yfo, zfo, cdata.orange.framenumber]];
end
% sort localization by time identified
sort_locs = sortrows(unsort_locs,3);
chunks = ceil(numel(unsort_locs(:,1))/mol_chunks); % we want a ceiling because we'd like to place all molecules in the data set

% initialize variables of interest
drift = [];
sub_locs_1 = sort_locs(1:mol_chunks,:);

% we're always grabbing i through i+1, so at final interation (chunks -1)
% we're going through (chunks-1)*mol_chunks : chunks*mol_chunks/end, which
% should be the total number of molecules in the dataset
for i = 1:chunks-1
    % Swap between a complete mol_chunk calculation or grab remaining
    try
        sub_locs_2 = sort_locs(i*mol_chunks+1:(i+1)*mol_chunks,:);
    catch
        sub_locs_2 = sort_locs(i*mol_chunks+1:end,:);
    end
    % impose previously detected drift on localization collection
    try
        sub_locs_2(:,1) = sub_locs_2(:,1) + sum(drift(:,1));
        sub_locs_2(:,2) = sub_locs_2(:,2) + sum(drift(:,2));
        sub_locs_2(:,3) = sub_locs_2(:,3) + sum(drift(:,3));
    catch
    end
    
    % ICP calculation
    d_max = 1; % Initial max distance should include all points
    % Reset cumulative drift variables
    x_drift = 0;
    y_drift = 0;
    z_drift = 0;
    % Iterative loop
    last_cost = -10000;
    j = 0;
    while true
        j = j+1;
        while true % iterative approach computes 'matches' of point clouds based on statistical information
            [I,distances] = knnsearch(sub_locs_1(:,1:3),sub_locs_2(:,1:3)); % nearest neighbor search
            % We want to always find the closest points in sub_locs_1 to
            % sub_locs_2. As we grow sub_locs_1, we'll have better sampling,
            % and sub_locs_2 will always stay about the same size
            
            included = 0;
            dx = 0;
            dy = 0;
           
            mu_distance = mean(distances(distances < d_max));
            sig_distance = std(distances(distances < d_max));
          
            if mu_distance < d_max % really good registration
                d_max = mu_distance + 3*sig_distance;
            elseif mu_distance < 3*d_max % good registration
                d_max = mu_distance + 2*sig_distance;
            elseif mu_distance < 6*d_max % ok registration
                d_max = mu_distance + sig_distance;
            else % terrible registration
                d_max = 1;
            end
            if abs(last_cost - mu_distance) < 0.00001 % Last cost will represent the previous average NN distance
                break; % If difference between iterations is sufficiently small, break the loop
            else
                last_cost = mu_distance;
            end
            
            
        end
        
        % Drift recording
        
        index = distances < d_max;
        data_0 = [sub_locs_1(I(index),1:3)].';
        
        data_1 = [sub_locs_2(index,1:3)].';
        
        [R, T] = dual_quaternion_transform(data_0,data_1,exp(-(distances(index)./(2*std(distances(index))))));
        
        x_drift = x_drift + T(1);
        y_drift = y_drift + T(2);
        z_drift = z_drift + T(3);
        
        sub_locs_2(:,1:3) = (sub_locs_2(:,1:3).' + T).';
        
        % While loop break conditions are convergence on correction vector
        % or exceeding maximum iterations
        if abs(T(1)) < 0.001 && abs(T(2)) < 0.001 && abs(T(3)) < 0.001
            break
        elseif j > max_iter
            break
        end
    end

    % Add corrected localizations to current 
    sub_locs_1 = [sub_locs_1;sub_locs_2];
    drift(i,:) = [x_drift,y_drift, z_drift, mean(sub_locs_2(:,4))];
end

tot_drift = [];
for i = 1:numel(drift(:,1))
    tot_drift(i,:) = [sum(drift(1:i,1)), sum(drift(1:i,2)), sum(drift(1:i,3)), drift(i,4)];
end

try % Try correcting red localizations if there are any
drift_correct_x = spline(tot_drift(:,4),gausssmooth(tot_drift(:,1),3,6),cdata.red.framenumber);
drift_correct_y = spline(tot_drift(:,4),gausssmooth(tot_drift(:,2),3,6),cdata.red.framenumber);
drift_correct_z = spline(tot_drift(:,4),gausssmooth(tot_drift(:,3),3,6),cdata.red.framenumber);
cdata_out.red.xf = xfr + drift_correct_x;
cdata_out.red.yf = xfr + drift_correct_y;
cdata_out.red.zf = xfr + drift_correct_z;
catch
end
try % try correcting orange localizations
drift_correct_x = spline(tot_drift(:,4),gausssmooth(tot_drift(:,1),3,6),cdata.orange.framenumber);
drift_correct_y = spline(tot_drift(:,4),gausssmooth(tot_drift(:,2),3,6),cdata.orange.framenumber);
drift_correct_z = spline(tot_drift(:,4),gausssmooth(tot_drift(:,3),3,6),cdata.orange.framenumber);
cdata_out.orange.xf = xfo + drift_correct_x;
cdata_out.orange.yf = yfo + drift_correct_y;
cdata_out.orange.zf = zfo + drift_correct_z;
catch
end