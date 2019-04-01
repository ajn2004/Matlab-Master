function [im1, um] = make_pics(fdname, varargin)
% This is a rednering function that will 
%code adapted from show_pics in da_support folder
load(fdname); 
numvar = length(varargin);
if numvar < 1
    sm = 1;
elseif numvar >= 1
   sm = varargin{1};
end
%% Figuring out pixel size
q = q*1000; % convert q to nm
xf_crlb = crlbs(:,1);
yf_crlb = crlbs(:,2);
xf_all = xf_fixed;
yf_all = yf_fixed;
zf_all = ncoords(:,3);

rad = mean((xf_crlb.*yf_crlb).^0.25)*q; % all distances in this are determined off of loc_unc
grid_size = 0.5*rad; % pixwidth in nm/pix
mag = q/grid_size; % Sam's 'expansion' factor, effectively zoom between pix and grid
um = 1000/grid_size; % number of pixels required to make a scale bar of 1 um
if numvar >= 2
    gride = varargin{2};
    if gride > 0
        grid_size = gride;
    end       
end
%% Intializing image grid
xmax = grid_size*ceil(max(xf_all)*q/grid_size); % Determine maximum x pixels
ymax = grid_size*ceil(max(yf_all)*q/grid_size); % Determine maximum y pixels
[Xgrid, Ygrid] = meshgrid(0:grid_size:xmax, 0:grid_size:ymax); % create x and y meshgrid
i1 = zeros(round(ymax/grid_size),round(xmax/grid_size),1); %intialize grayscale image variable


%% Populating image grid
% Color will be selected from a pallet and added to an rgb image as
% described above
for i = 1:numel(xf_all)
    x_ind = find(Xgrid(1,:) > xf_all(i)*q, 1, 'first') - 1;
    y_ind = find(Ygrid(:,1) > yf_all(i)*q, 1, 'first') - 1;
    
    % Set image RGB
    i1(y_ind,x_ind,1) = i1(y_ind,x_ind,1) + 1;

end

% Blur the image
[x,y] = meshgrid(-50:50,-50:50);
G = exp(-(x.^2 + y.^2)./(0.5*(sm*rad)^2));
gn = G./(sum(G(:)));

for i = 1:1
    im1(:,:,i) = conv2(i1(:,:,i),gn,'same');
end

im1(:,:,1) = im1(:,:,1)./max(max(im1(:,:,1)));
clear sim1 G gn x y drifts fpath framenum_all i iln iloc imgf llv N N_crlb nfiles off_all
clear off_crlb pixtpho pixw q sigx_all sigx_crlb sigy_all sigy_crlb simgf
clear total_mol total_molecues xdrift xf_all xf_crlb xf_fixed xf_part ydrift yf_all yf_crlb
clear yf_fixed yf_max zf_nm zf_nm_temp
save(['render\',fdname(1:end-4),'_render.mat']);