%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background-Away
%
% This will be the basis for a new method of background estimation and
% subtraction based on either wavelet transformation and removal of low
% frequency 3D information or through spline fitting of background
%
%
% AJN 5/26/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars; close all; clc;
[fname, fpath] = uigetfile('*tif');
cd(fpath);
i1= readtiff(fname);
figure('Units','Normalized','Outerposition',[0,0,0.5,1])
i1s = i1;
[m,n,p] = size(i1s);
count = 1;
% for sp = 0:100
sp = 10;
    i1s = i1;
%     clear i3
    i2 = fft(i1s,p,3);
    clear i1s
    i2(:,:,1:sp) = i2(:,:,1:sp)*0;
    i2(:,:,end-sp+1:end) = i2(:,:,end-sp+1:end)*0;
    i3 = ifft(i2,p,3);
    clear i2
    imagesc(mean(abs(i3),3))
    colormap('gray')
    title(['sp = ', num2str(sp)]);
    drawnow
    M(count) = getframe(gcf);
    count = count+1;
% end