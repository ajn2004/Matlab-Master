grid = 30;
[Rx, Ry, Rz] = rot_mat(pi/2+1);
% for l = 0:360
[R2x, R2y,R2z] = rot_mat(pi/2 + 2.01);
sm = 5;
t = (R2y*Rx*[xs,ys,zs].').';
xr = t(:,1);
yr = t(:,2);
zr = t(:,3);
xxs = xr - min(xr);
yys = yr - min(yr);

%% Intializing image grid
xmax = grid_size*ceil(max(xxs)/grid_size); % Determine maximum x pixels
ymax = grid_size*ceil(max(yys)/grid_size); % Determine maximum y pixels
[Xgrid, Ygrid] = meshgrid(0:grid_size:xmax, 0:grid_size:ymax); % create x and y meshgrid
i1 = zeros(round(ymax/grid_size),round(xmax/grid_size),1); %intialize grayscale image variable


%% Populating image grid
% Color will be selected from a pallet and added to an rgb image as
% described above
for i = 1:numel(xxs)
    x_ind = find(Xgrid(1,:) > xxs(i), 1, 'first') - 1;
    y_ind = find(Ygrid(:,1) > yys(i), 1, 'first') - 1;
    
    % Set image RGB
    i1(y_ind,x_ind,1) = i1(y_ind,x_ind,1) + 1;

end

% Blur the image
[x,y] = meshgrid(-50:50,-50:50);
G = exp(-(x.^2 + y.^2)./(0.5*(sm*rad)^2));
gn = G./(sum(G(:)));
clear im1
for i = 1:1
    im1(:,:,i) = conv2(i1(:,:,i),gn,'same');
end
imagesc(im1)
colormap('gray')
% M(l+1) = getframe(gcf);
% end