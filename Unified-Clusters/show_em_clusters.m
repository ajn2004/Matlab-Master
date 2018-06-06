% Cluster_tolerance
clearvars;
close all;
clc;

min_num = 40;
min_area = 0;
min_perim = 0;


[fname, fpath] = uigetfile('*meas.mat');
load([fpath,fname]);
count = 1;
for i = 1:numel(c)
    if (c(i).mem_num >= min_num && c(i).area >= min_area && c(i).perimeter >= min_perim)
        scatter(c(i).sub_data(:,1),c(i).sub_data(:,2),'filled')
        hold on
        plot(c(i).orbit(:,1),c(i).orbit(:,2))
        hold off
        title(['Cluster ', num2str(i), ' of ', num2str(c(i).mem_num),' localizations']);
        xlabel('Position in um');
        ylabel('Position in um');
        axis image
%         pause(0.5)
       drawnow
       area(count) = c(i).area;
       perim(count) = c(i).perimeter;
       mems(count) = c(i).mem_num;
       M(count) = getframe(gcf);
            count = count +1;
        
    end
    
end
% writetiff(i1,[fname(1:end-4),'_clusters.tif']);
movie2gif(M,'cluster_animation.gif','DelayTime', 0.57);