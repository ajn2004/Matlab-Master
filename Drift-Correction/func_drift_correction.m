function corrections = func_drift_correction(data, pix_size, chunk_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function Drift Correction 
%
%  A function that analyzes XY-T localization data and return
%  AJN 7-21-16
%  Assume data = [ncoords(:,1),ncoords(:,2),framenum_all(:)]
%
%
% 

p = pix_size/1000;
gauss_std = 1.75; % width of gaussian to convolve with in pixels

% Data welding
framenum_all = data(:,3);
xf_all = data(:,1);
yf_all = data(:,2);
xf_in = xf_all;
yf_in = yf_all;
max_frames = framenum_all(end);
corrections = [xf_all,xf_all]*0;
% max_frames = max(framenum_all);
chunks = floor((max_frames -1)/chunk_size)+1;
[X, Y] = meshgrid(-10:10,-10:10);
gauss = exp(-2*((X).^2 + (Y).^2)./(gauss_std*2)^2);
% Build the images for subsequent correlations
for j = 1: chunks 
    % build xf and yf part for density plot rendering
    frame1 = (j-1)*chunk_size + 1;
    frame2 = j*chunk_size+1;
    ind = framenum_all >= frame1 & framenum_all <= frame2;
    xf_part = xf_in(ind);
    yf_part = yf_in(ind);
    
    % determine maximum plotted position in microns
    max_x = ceil(max(xf_in)/p)*p+p;
    max_y = ceil(max(yf_in)/p)*p+p;
    
    [Xgrid, Ygrid] = meshgrid(0:p: max_x,0:p: max_y);
    dens = zeros(size(Xgrid));
    
    
    [m, n] = size(Xgrid);
    for i = 1:numel(xf_part)
        x_ind = find(Xgrid(1,:) > xf_part(i), 1, 'first') - 1;
        y_ind = find(Ygrid(:,1) > yf_part(i), 1, 'first') - 1;
        dens(y_ind,x_ind) = dens(y_ind,x_ind) + 1;
    end

    
    dens = conv2(dens,gauss,'same');
    cmd_str = ['im',num2str(j),'=dens;'];
    eval(cmd_str);
%     cmd_str = ['imagesc(im',num2str(j),');'];
%     eval(cmd_str);
%     pause(0.2);
end
[m,n] = size(im1);
% perform cross correlations and find shift
drifts = [0,0];
for i = 1:chunks-1 % start on frame i and correlate with frame i+1
    cmd_str = ['the_coords = xcorr2(im',num2str(1), ', im',num2str(i+1),');'];
    eval(cmd_str);
%     sub_coords = the_coords(m-10:m+10,n-10:n+10);
    imgf = imgaussfilt(the_coords,2);
    simgf = imgf(m-5:m+5,n-5:n+5);
    [row, col] = find(imgf == max(max(simgf)));
    drifts = [drifts; pix_size*(col(1) - n), pix_size*(row(1) - m)];
end

% correct local drift to absolute drift
for i = 2:numel(drifts(:,1))
    drifts(i,:) = drifts(i-1,:) + drifts(i,:);
end

% determining functional form of drift in x and y
%in this section we will loop over a 9th order polynomial 
n = numel(drifts(:,1));
x = (0.5:1:chunks-0.5)*chunk_size;
x = x.';
% Calculate fit through loraxian transform
% Calculate fit through loraxian transform
% yout = func_lorax(x,drifts,framenum_all,100,0.01, 18, 10000);
% yout(1:2,:) = 0*yout(1:2,:);

driftx  = gausssmooth(drifts(:,1), 1.5, 4);
drifty  = gausssmooth(drifts(:,2), 1.5, 4);
sx = spline(x,driftx,framenum_all);
sy = spline(x,drifty,framenum_all);
yout = [sx,sy];
% display fit
figure
subplot(2,2,1);
plot(x,drifts(:,1),'.b');
hold on
plot(framenum_all,yout(:,1));
hold off
xlabel('frame chunks'); 
ylabel('drift in nm'); 
title('X drift');

subplot(2,2,3);
plot(x,drifts(:,2),'.b');
hold on
plot(framenum_all,yout(:,2));
xlabel('frame chunks'); 
ylabel('drift in nm'); 
title('Y drift');

subplot(2,2,[2,4]);
plot(drifts(:,1),drifts(:,2));
xlabel('xdrift');
ylabel('ydrift');
title('Total Drift');
axis equal

% correct data
xdrift = yout(:,1)/(p);
ydrift = yout(:,2)/(p);
corrections = [xdrift, ydrift];
xf_fixed = xf_all+xdrift;
yf_fixed = yf_all+ydrift;

plot(xf_all,yf_all,'.r');
hold on
plot(xf_fixed,yf_fixed,'.b');
hold off
legend('Original','Fixed');
