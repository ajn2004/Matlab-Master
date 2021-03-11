function cdata = model_drift_correction(cdata, align_color, frames_in_chunk )

% Build initial coordinate set and overlay via neural net transform
% xfo and xfr were initially to indicate orange from red, but given the
% current analysis is color indiscriminante this is a vestige from
% development that will be kept as xfr: reference distribution xfo: off
% color distribution
[xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);

if strcmpi(align_color,'orange') % If we specify orange, make reference distribution be orange molecules
    [xfr, yfr] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
    zfr = -cdata.orange.zf;
    xfo = cdata.red.xf;
    yfo = cdata.red.yf;
    zfo = -cdata.red.zf;
    frame_o = cdata.red.framenumber;
    frame_r = cdata.orange.framenumber;
    xfo = xfo-mean(xfr);
    zfo = zfo-mean(zfr);
    yfo = yfo-mean(yfr);
    xfr = xfr-mean(xfr);
    zfr = zfr-mean(zfr);
    yfr = yfr-mean(yfr);
    next_color = 'red';
else
    [xfo, yfo] = make_nn_channel_transform(cdata.orange.xf,cdata.orange.yf);
    zfo = -cdata.orange.zf;
    xfr = cdata.red.xf;
    yfr = cdata.red.yf;
    zfr = -cdata.red.zf;
    xfo = xfo-mean(xfr);
    zfo = zfo-mean(zfr);
    yfo = yfo-mean(yfr);
    xfr = xfr-mean(xfr);
    zfr = zfr-mean(zfr);
    yfr = yfr-mean(yfr);
    frame_r = cdata.red.framenumber;
    frame_o = cdata.orange.framenumber;
    next_color = 'orange';
end

% Make data arrays
r_data = [xfr, yfr, zfr];
o_data = [xfo, yfo, zfo];

% Build t0 model by segmenting out early data
index = frame_r >= 0 & frame_r < frames_in_chunk;
data0 = r_data(index,:);
epsb = 0.2;
epsn = 0.006;
amax = 50;
lamb = 100;
alpha = 0.5;
decr = 0.995;
vangle = [00, 90];
params = [epsb, epsn, amax, lamb, alpha, decr];
[w, T, A] = func_gng(data0, params);
w0 = w;
% Reduce w0 down to a decreased number of nodes
% Every decrease, replace existing nodes using a k-means / radius
radius = 0.1; % micron radius for consideration
seen = w0(:,1)*0; % index variable to record whether a node has been removed
w1 = [];
for i = 1:numel(w0(:,1)) % loop over every node
    if ~seen(i) % Only look at new molecules
        seen(i) = 1; % record seeing this molecule
        w_dist = w0(:,1)*0;
        d_dist = data0(:,1)*0;
        for k = 1:3
            w_dist = (w0(:,k) - w0(i,k)).^2 + w_dist;
            d_dist = (data0(:,k) - w0(i,k)).^2 + d_dist;
        end
        % look to see if any nodes are too close
        w_count = sum(w_dist.^0.5 < radius)-1; % We have to correct for 1 because ith data set is included
        d_count = sum(d_dist.^0.5 < radius);
        if d_count > 5
            if w_count >= 1 % If there are multiple overlaping points, merge them into one
                ind = w_dist.^0.5 < radius & seen == 0; % grab all points overlapping
                seen(ind) = 1;
                w1 = [w1;[mean(w0(ind,1)),mean(w0(ind,2)),mean(w0(ind,3))]];
            else
                w1 = [ w1;w0(i,:)];
            end
        end
    else
    end
end

for l = 1:10
    for i = 1:numel(w1(:,1))
        d_dist = data0(:,1)*0;
        for k = 1:3
            d_dist = (data0(:,k) - w1(i,k)).^2 + d_dist;
        end
        index = d_dist.^0.5 < radius;
        w1(i,:) = mean(data0(index,:),1);
    end
end


for i = 1:100 % perform individual step  
        IDx = knnsearch(w1,data0);
        for i = 1:numel(w1(:,1))
            index = IDx == i; % Grab all points whose NN is node i
%             Determine error due to data position
            if sum(index) > 1
                       
            for k= 1:3
                dx(k) = mean(data0(index,k)) - w1(i,k); % This is the correction value based on the data
                w1(i,k) = w1(i,k) +dx(k);
            end
            end
    
        end
        
end
keep_index = find(w1(:,1) == w1(:,1));
w0 = w1(keep_index,:);
% Upd cdata.red.framenumber*0ate neural points to new distribution

node_association =[];
onode = [];
odata = [xfo,yfo,zfo];
data = [xfr, yfr, zfr];

clear w_all
w_all(:,:,1,1) = w0; % last index = 1 is previous frames model adapted to current time chunk
w_all(:,:,1,2) = w0; % last index = 2 is current time chunk's adapted model
% We are going to find the average drift correction between frames,
% Every frame will then need 2 different nodal sets, the previous frames
% node and the current frame's node. 
for i = 2:10
    r_index = frame_r >= (i-1)*frames_in_chunk & frame_r < i*frames_in_chunk;
    o_index = frame_o >= (i-1)*frames_in_chunk & frame_o < i*frames_in_chunk;
    data1 = data(r_index,:);
    data2 = odata(o_index,:);
    % Data 1 and data 2 are time separated out localization sets
    
    % Change the average position of the nodes
    % Grab previous frames adapted model for average comparison
    w_all(:,:,i,1) = update_average_neural_gas(w_all(:,:,i-1,2), data1, 0.2);
    node_association = [node_association;knnsearch(w_all(:,:,i,1),data1,"K",3)];
    onode = [onode;knnsearch(w_all(:,:,i,1),data2,"K",3)];
    w_all(:,:,i,2) = update_local_neural_gas(w_all(:,:,i,1), data1, 0.2);

end

% Stitch together the average cell drift
drift = [0, 0, 0];
dw = [;]
for i = 2:10
    dw = mean(w_all(:,:,i,1) - w_all(:,:,i-1,2)); % Take difference from identical drift models
    drift(i,:) = drift(i-1,:) + dw;
end
% Correct the drift
data_c = data;
data_o = odata;
for i = 1:3
    red_model = spline(500:1000:9500, drift(:,i),frame_r);
    orange_model = spline(500:1000:9500, drift(:,i),frame_o);
    data_c(:,i) = data(:,i) - red_model;
    data_o(:,i) = (odata(:,i) -orange_model);
    
end
data_c(:,3) = data_c(:,3)*2.5;
data_o(:,3) = data_o(:,3)*2.5;

cdata.(align_color).xf = data_c(:,1);
cdata.(align_color).yf = data_c(:,2);
cdata.(align_color).zf = data_c(:,3);
cdata.(next_color).xf = data_o(:,1);
cdata.(next_color).yf = data_o(:,2);
cdata.(next_color).zf = data_o(:,3);
end