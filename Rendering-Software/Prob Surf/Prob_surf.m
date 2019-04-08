% Build 3D probability density map
clearvars;
close all;
clc;

% This algorithm will perform by building a 3D histogram of cube size 'bin'
% From there it will perform a 'fluorescence' ray tracing onto a pixel map
% defines as m x n

[m,n] = hd_stats(720);
m = 245;
n = 268;
bin = 70; % bin size in nm



% Choose Image
fname = 'hek6_r2_dz20_dast_tol_dc_100nm_traj.mat';
load(fname);
r = str2num(fname(strfind(fname,'_r')+2));
zf = func_shift_correct(ncoords(:,3)*q,framenumber,r);
zf = zf(:)-mean(zf);
try
    xf = q*xf_fixed;
    yf = q*yf_fixed;
    
catch lsterr
    xf = q*ncoords(:,1);
    yf = q*ncoords(:,2);
end
xf = xf - mean(xf);
yf = yf - mean(yf);

for p = 0
    [rx,ry,rz] = rot_mat(deg2rad(10*p));
    coords = [xf,yf,zf];
    rcoords = (ry*coords.').';
i1 = func_3D_hist(rcoords,bin);
% clearvars -except i1 m n fname
[m0,n0,o0] = size(i1);
% for i = 1:numel(i1(1,1,:))
%     imagesc(i1(:,:,i));
%     axis image
%     drawnow
% end
imagesc(sum(i1,3))
drawnow

% At this point a 3D histogram has been created that contains the
% localization population data Now to perform the Ray-Tracing calculation

% Building final image grid
i2 = zeros(m,n);
mmag = (m)/(m0);
nmag = (n)/(n0);
mag = max([nmag,mmag]);
L = -100;
lp = -L/mag;

zsp = (1:o0) + lp+100;
zs = 1:o0;
% X is being defined as the column-coordinate of the output image
% Y is defined as the row-coordinate
z = L; % Arbitrary value as the 'mag' scales this


for i  = 1:m % calculate the intensity at each image pixel
    tic
    y = i-m/2; % Row position of pixel being considered
    my = atan(y/z); % Determine the 'y-slope' of the line to intersect object pixels
    for j = 1:n
        % Define pixel of interest's initial coord
        x = j-n/2;
        mx = atan(x/z); % Determine the 'x-slope' of the line to intersect object pixels
        xs = round(mx*zsp + n0/2)+1;
        ys = round(my*zsp + m0/2)+1;
        
        for k = 1:numel(zsp)
            try
                i2(i,j) = i2(i,j) + i1(ys(k),xs(k),zs(k));
            catch
            end
        end
%         ind = find(abs(islope(:,:,:,1)-my) <= dy & abs(islope(:,:,:,2) - mx) <=dx);
%         i2(i,j) = sum(sum(i1(ind)));
    end
    t(i) = toc;
    ajn_wait(t,i,m)
end
imagesc(i2);
% M(p) = getframe(gcf);
end