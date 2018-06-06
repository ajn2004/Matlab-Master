function func_time_series_2(ifin, peakresponse, psize, l,sn)
close all


flag = 0;
%% Manual Select of measuring regions
if exist(['measured_results',num2str(l-2),'.mat'])
    load(['measured_results',num2str(l-2),'.mat'],'points');
    flag = 1;
else
    points = [];
end

if strcmp(sn,'y') || strcmp(sn,'Y')
    load(['measured_results',num2str(l-1),'.mat']);
else
    figure('units','normalized','outerposition',[0 0 0.5 1]);
    w = 1;
    im2show = var(ifin,1,3);
    im2 = (im2show >0).*im2show;
    imagesc(im2, [mean(im2(:)), max(im2(:))/2]);
    axis image
    colormap('jet')
    colorbar
    if numel(points) > 0
        hold on
        plot(points(:,1),points(:,2),'x','MarkerSize',4)
        
        hold off
        
    end
    pl = gca;
    clear points
    hold on
    while true
        clearvars cents;
        % selects point clicked in plot
        title('Select a new center, press enter to quit')
        bc = waitforbuttonpress;
        if bc == 1
            break
        else
            cents = get(pl,'currentpoint');
            % assigns selected points to array
            points(w,1) = round(cents(1,1)); % save value of center in um
            points(w,2) = round(cents(1,2)); % save value of center in um
            
            
            makebox(points(w,:),psize);
            w = w+1;
        end
    end
    hold off
end
if numel(points) > 0
    [fluor, xs] = func_measure_regions(ifin./peakresponse, points, psize);
    save(['measured_results', num2str(l-1),'.mat'],'xs','points','fluor');
end