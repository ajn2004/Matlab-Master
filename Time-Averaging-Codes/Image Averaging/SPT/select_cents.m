function [points] = select_cents(i1,col, psize)
points = [];
% A script to select points
% Initialize the figure
figure('units','normalized','outerposition',[0 0 1 1]); 
w = 1;
imagesc(i1);
colormap(col)
colorbar
pl = gca;
clear points
hold on
% start selecting points
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
        %         plot(cents(:,1),cents(:,2),'.r','MarkerSize',35)
        
        makebox(points(w,:),psize);
        w = w+1;
    end
end
end