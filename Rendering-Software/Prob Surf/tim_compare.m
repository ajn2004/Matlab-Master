%% Speed Test
close all;
clearvars;
clc;

cd('C:\Users\AJN Lab\Documents\GitHub\Matlab-Master\Rendering-Software\Prob Surf');
load('C:\Users\AJN Lab\Dropbox\Data\4-2-19 hek-3d-trial\Analysis\toleranced4thru8\DC\hek5_r2_dz20_dast_tol_dc.mat')
xf = xf_fixed*q;
yf = yf_fixed*q;
zf = ncoords(:,3)*q;
zf = func_shift_correct(ncoords(:,3)*q,framenumber,2).';
xf = xf - mean(xf);
yf = yf - mean(yf);
zf = zf - mean(zf);
[rx,ry,rz] = rot_mat(deg2rad(0));
rots = rx*[xf,yf,zf].';
% sizes = [500, 400, 300, 200, 100, 90, 80, 70, 60, 50, 40, 30, 25]; 
% for kk = 1:numel(sizes)
for i = 1:100
    tic
% disp('GOING!')
% try
    i1 = func_3D_dens(single(rots).',30,0.06);
% catch lsterr
% end
    gpu_t(i) = toc
% toc
% %     
end

% [m,n,o] = size(i1);
% pixes(kk) = m*n*o;
% gputim(kk) = mean(gpu_t);
% for tm = 1:100
%     tic
%     i2 = func_3D_hist(rots.',sizes(kk),0.06);
%     cpu_t(tm) = toc;
%     ajn_wait(cpu_t,kk,100);
% end
% cputim(kk) = mean(cpu_t);
% end

%     for i = 1:o
%         imagesc(i1(:,:,i));
imagesc(mean(i1,3));
%     imagesc(i1(:,:,2));
        axis image
        drawnow
%         M(i) = getframe(gcf);
%     end
%     for i = numel(M):-1:1
%         imagesc(i1(:,:,i));
%         axis image
%         drawnow
%         M(numel(M)+1) = getframe(gcf);
%     end
%     movie2gif(M,'diffraction_roll_image.gif','Delaytime', 0.05,'LoopCount',Inf);
%     plot3(rots(1,:),rots(2,:),rots(3,:),'.')
%     xlim([-10 10])
%     ylim([-10 10])
%     zlim([-10 10])
%     drawnow
% end