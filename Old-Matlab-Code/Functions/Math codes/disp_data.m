function [avef, num] = disp_data(gs)

tinfo = dir('*results*');
% figure('units','Normalized','OuterPosition',[0,0,1,1]);
% hold on
totf = [];
num = 0;
for i = 1:numel(tinfo)
    load(tinfo(i).name);
    totf = [totf,fluor(:)]; 
    if exist('xs')
    num = num + numel(xs(:,1));
    end
end
avef = mean(totf,2);
% subplot(1,2,1);
% errorbar(1:numel(avef),avef,std(totf,1,2)./numel(totf(1,:))^0.5)
% plot(avef);
% hold off
% title('dF/dF_m (error bars are SOM');
% xlabel('Frame number');
% ylabel('Average normalized differential fluorescent response');
% 
plot(gausssmooth(avef,gs,10));
title('Gaussian Smoothed Delta F over F')
xlabel('Frame')
ylabel('Average normalized differential fluorescent response');
hold on
for i = 50:10:100
    plot([i, i],[-1,1],'r');
end
% ylim([-0.001,.006]);
% xlim([45,100]);
legend(['1/e^2 smoothing radius = ',num2str(gs),' frames'],'Stimulus');
hold off
% subplot(1,2,2);
% df = diff(gausssmooth(avef,gs,10));
% plot(df);
% title('Derivative of Gaussian Smoothed Delta F over F');
% xlabel('Frame');
% ylabel('Change in fluorescent signal');
% ylim([-0.01,0.01]);
end