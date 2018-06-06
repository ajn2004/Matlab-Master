function peakresponse = getpeakresponse(fname, psize)
%GETPEAKRESPONSE get the average peak response from a data set exposed to
%ammonium chloride
%   A = GETPEAKRESPONSE(FILENAME) will return a value which is the average
%   peak phluorin response which can is delta F_NH4CL
%
    %% Load Baseline image
    i1 = readtiff(fname);
    iprod = rollingball(i1,4,4);
    s1 = sum(sum(i1)); % sum all frames and store the final result
    ind = find(s1 == max(s1));
    sm1 = gausssmooth(s1,5,20);
    plot(sm1(25:end-25));
    title('Select a peak response point')
    [x1, y1] = ginput(1);
    imagesc(iprod(:,:,round(x1)));
    pl = gca;
    clear points
    hold on
    w=1;
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
    hold off

    [fluor, xs] = func_measure_regions(iprod, points, psize);
    smfluor = gausssmooth(fluor,5,20);
    plot(smfluor(20:end-20));
    title('Select 3 points in this order [base_end, peak_start, peak_end]');
    [x,y] = ginput(3);
    if x(2) < x(3)
        xlow = round(x(2));
        xhigh = round(x(3));
    else
        xlow = round(x(3));
        xhigh = round(x(2));
    end
    
    peakresponse = mean(fluor(1,1,xlow:xhigh)-mean(fluor(1,1,1:x(1))));
%     peakr = 50;
    
%     peakresponse = i1(:,:,peakr);
    
end