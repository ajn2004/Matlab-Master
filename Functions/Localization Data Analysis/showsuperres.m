function im1 = showsuperres(xf_all,yf_all,varargin)
num = numel(xf_all);
if numel(varargin) >= 1
    q = varargin{1};
    xf_in = xf_all*q;
    yf_in = yf_all*q;
    p = (max(xf_in) - min(xf_in))/ numel(xf_all)^0.5;
    gauss_std = 2;
    if numel(varargin) >= 2
        p = varargin{2};
    end
    if numel(varargin) >= 3
        gauss_std = varargin{3};
    end
    if numel(varargin) == 5
        max_x = varargin{4};
        max_y = varargin{5};
    end
else
    xf_in = xf_all;
    yf_in = yf_all;
end

% max_x = ceil(max(xf_in)/p)*p+p;
max_x = 80;
% max_y = ceil(max(yf_in)/p)*p+p;
max_y = 80;
[Xgrid, Ygrid] = meshgrid(0:p: max_x,0:p: max_y);
dens = zeros(size(Xgrid));
[m, n] = size(Xgrid);
for i = 1:num
    x_ind = find(Xgrid(1,:) > xf_in(i), 1, 'first') - 1;
    y_ind = find(Ygrid(:,1) > yf_in(i), 1, 'first') - 1;
    dens(y_ind,x_ind) = dens(y_ind,x_ind) + 1;
end

% imagesc(dens)
[X, Y] = meshgrid(-10:10,-10:10);
gauss = exp(-2*((X).^2 + (Y).^2)./(gauss_std*2)^2);

im1 = conv2(dens,gauss,'same');
% w = 1;
% imagesc(im1);
% drawnow
end