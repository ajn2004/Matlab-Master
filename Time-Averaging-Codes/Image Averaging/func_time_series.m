function [fluor, points, xs] = func_time_series(ifin, peak, ave2,peakresponse, psize, baseline)
close all



flag = 0;
%% Manual Select of measuring regions
if exist('measured_results.mat')
    ans = input('A set of measured values has been found, would you like to use them? (y/n)','s');
    if strcmp(ans,'Y') || strcmp(ans,'y')
        flag =1;
    end
end
if flag ==0
    [points] = func_buton_id(ifin, baseline);
    points = round(points);
    for j = 1:4
    for i = 1:numel(points(:,1))
        [row, col] = find(sum(ifin(points(i,2)-psize :points(i,2)+psize ,points(i,1)-psize :points(i,1)+psize ,:),3).^2 == max(max(sum(ifin(points(i,2)-psize :points(i,2)+psize ,points(i,1)-psize :points(i,1)+psize ,:),3).^2)));
        points(i,:) = [points(i,1) + (col(1) -3),points(i,2) + (row(1)-3)];
    end
    end
    figure('units','normalized','outerposition',[0 0 1 1]);
    w = 1;
    imagesc((sum(ifin,3)>0).*sum(ifin,3).^2);
    
    colormap('hsv')
    colorbar
    hold on
    plot(points(:,1),points(:,2),'x')
    hold off
    axis image
    waitforbuttonpress

%     pl = gca;
%     
%     hold on
%     while true
%         clearvars points;
%         % selects point clicked in plot
%         title('Select a new center, press enter to quit')
%         bc = waitforbuttonpress;
%         if bc == 1
%             break
%         else
%             points = get(pl,'currentpoint');
%             % assigns selected points to array
%             cents(w,1) = round(points(1,1)); % save value of center in um
%             cents(w,2) = round(points(1,2)); % save value of center in um
%             %         plot(cents(:,1),cents(:,2),'.r','MarkerSize',35)
%             
%             makebox(cents(w,:),psize);
%             w = w+1;
%         end
%     end
%     hold off
else
    load('measured_results.mat');
end
% Growing Neural Gas
% [cents] = func_grow_neural_gas(ifin(:,:,peak));
[fluor, xs] = func_measure_regions(ifin./peakresponse, points, psize);
save('measured_results.mat','points','fluor');
plot(100*fluor(:))
title('Averaged Time Series of dF/max(F_N_H_4_C_l))');
xlabel('Frame Number');
ylabel('Normalized Fluorescence(%)');