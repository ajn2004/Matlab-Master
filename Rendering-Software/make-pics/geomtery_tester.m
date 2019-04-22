%% Speed Test
close all;
clearvars;
clc;

bin = [1080,1920];
pnts = 4000;
% cd('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Rendering-Software\Prob Surf');
load('C:\Users\AJN Lab\Dropbox\Data\4-2-19 hek-3d-trial\Analysis\toleranced4thru8\DC\hek5_r2_dz20_dast_tol_dc.mat')
xf = xf_fixed*q;
yf = yf_fixed*q;
zf = ncoords(:,3)*q;
zf = func_shift_correct(ncoords(:,3)*q,framenumber,2).';
%% Generate Simulation Data
% R = 1; % Radius in microns
% for i = 1:pnts
%     theta = rand*2*pi;
%     phi = rand*pi;
%     zf(i,1) = rand;
%     xf(i,1) = R*cos(zf(i)*2*pi/0.3);
%     yf(i,1) = R*sin(zf(i)*2*pi/0.3);
%     
% %     xf(i,1) = rand*R;
% %     yf(i,1) = rand*R;
% %     zf(i,1) = rand*R;
%     
% end
xf = xf - mean(xf);
yf = yf - mean(yf);
zf = zf - mean(zf);
for i = 1:45
[rx,ry,rz] = rot_mat(deg2rad(i));
rots = rx*[xf,yf,zf].';
tic
    [i1] = make_pics_gpu([rots].',bin,1);
toc
[m,n,o] = size(i1);

imagesc(i1)
axis image
drawnow
M(i) = getframe(gcf);
end
for i = 1:1
[rx,ry,rz] = rot_mat(deg2rad(1));
rots = ry*[rots];
tic
    [i1] = make_pics_gpu([rots].',bin);
toc
[m,n,o] = size(i1);

imagesc(i1)
axis image
drawnow
M(numel(M)+1) = getframe(gcf);
end

xr = rots(1,:);
yr = rots(2,:);
zr = rots(3,:) - min(rots(3,:)) + 1;



