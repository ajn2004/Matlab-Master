% Build 3D probability density map
clearvars;
close all;
clc;

% This algorithm will perform by building a 3D histogram of cube size 'bin'
% From there it will perform a 'fluorescence' ray tracing onto a pixel map
% defines as m x n

[m,n] = hd_stats(720);
m = 100;
n = 100;
bin = 70; % bin size in nm



% Choose Image
fname = 'hek6_r2_dz20_dast_tol_dc_100nm_traj.mat';
load(fname);
r = str2num(fname(strfind(fname,'_r')+2));
zf = func_shift_correct(ncoords(:,3)*q,framenumber,r);
zf = zf(:);
try
    xf = q*xf_fixed;
    yf = q*yf_fixed;
    
catch lsterr
    xf = q*ncoords(:,1);
    yf = q*ncoords(:,2);
end

i1 = func_3D_hist([xf,yf,zf],bin);
clearvars -except i1 m n fname
[m0,n0,o0] = size(i1);
for i = 1:numel(i1(1,1,:))
    imagesc(i1(:,:,i));
    axis image
    drawnow
end
% At this point a 3D histogram has been created that contains the
% localization population data Now to perform the Ray-Tracing calculation

% Building final image grid
i2 = zeros(m,n);
mmag = m/m0;
nmag = n/n0;
mag = min([nmag,mmag]);
L = 100;
lp = L/mag;
[mgrid,zmgr] = meshgrid(1:m0,1:o0);
[ngrid,zngr] = meshgrid(1:n0,1:o0);
mg = mgrid - m0/2;
ng = ngrid - n0/2;
zmg = zmgr + lp - 1;
zng = zngr + lp - 1;
mms = atan(mg./zmg);
nms = atan(ng./zng);
islope(:,:,:,1) = i1*0;
islope(:,:,:,2) = i1*0;
% Build a 'slope' object for subsequent selection
for i = 1:m0
    for j = 1:n0
        for k = 1:o0
            islope(i,j,k,1) = mms(k,i);
            islope(i,j,k,2) = nms(k,j);
        end
    end
end
% X is being defined as the column-coordinate of the output image
% Y is defined as the row-coordinate
z = L; % Arbitrary value as the 'mag' scales this
dx = atan(0.5/z);
dy = atan(0.5/z);

for i  = 1:m % calculate the intensity at each image pixel
    tic
    y = i-m/2; % Row position of pixel being considered
    my = atan(y/z); % Determine the 'y-slope' of the line to intersect object pixels
    for j = 1:n
        % Define pixel of interest's initial coord
        x = j-n/2;
        mx = atan(x/z); % Determine the 'x-slope' of the line to intersect object pixels
        ind = find(abs(islope(:,:,:,1)-my) <= dy & abs(islope(:,:,:,2) - mx) <=dx);
        i2(i,j) = sum(sum(i1(ind)));
    end
    t(i) = toc;
    ajn_wait(t,i,m)
end
imagesc(i2);