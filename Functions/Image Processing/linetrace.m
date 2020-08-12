function [path,brightness] = linetrace(i1,varargin)
% This is a function to measure a line trace across image i1
% The option to prove a path in the form of [X,Y] 2-column vector is
% supported
if numel(varargin) == 0
    y = [];
    x = [];
    imagesc(i1);
    hold on
    
    pl = gca;
    clear points
    hold on
    w = 1;
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
            plot(points(:,1),points(:,2),'r');
            w = w+1;
        end
    end
    hold off
    path = [];
    for i = 2:numel(points(:,1))
        % find indices between points using polyfit command
        if points(i-1,1) <= points(i,1)
            xt = points(i-1,1):points(i,1);
        else
            xt = points(i,1):points(i-1,1);
        end
        if points(i-1,2) <= points(i,2)
            yt = points(i-1,2):points(i,2);
        else
            yt = points(i,2):points(i-1,2);
        end
        if numel(xt) >= numel(yt)
            a = polyfit(points(i-1:i,1),points(i-1:i,2),1);
            x = [x,xt];
            y = [y,round(a(1)*xt + a(2))];
        else
            a = polyfit(points(i-1:i,2),points(i-1:i,1),1);
            y = [y,yt];
            x = [x,round(a(1)*yt + a(2))];
        end
        
    end
    
    path = [x.',y.'];
elseif numel(varargin) == 1
    path = varargin{1};
end
for i = 1:numel(path(:,1))
    brightness(i) = i1(path(i,2),path(i,1));
end