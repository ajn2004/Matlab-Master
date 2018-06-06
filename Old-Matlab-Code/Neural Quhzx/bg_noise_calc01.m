function bkgn=bg_noise_calc01(imagefile,pix_to_pho,n_start,base_name)
% from Sam Hess lab, not written by AJN
    i1=imread(imagefile,n_start);
    h = figure('WindowStyle','modal','NumberTitle','off','Name',...
                'Select Region For Background Noise:',...
                'Pointer','cross','Units','Normalized','OuterPosition',[0 0 1 1]);
    set(h,'Units','Pixels');
    clims=[min(i1(:)) max(i1(:))];
    imagesc(i1,'HitTest','off',clims);
    title(base_name(1:end-4),'interpreter','none');
    axis image
    colormap gray
    drawnow
    try   %gets coordinates of user selected ROI
        waitforbuttonpress 
    catch
        return
    end
    p1 = get(gca,'CurrentPoint');
    rbbox;
    p2 = get(gca,'CurrentPoint');
    close(h);
    pause(.1); %makes sure figure closes

    p1 = double(uint16(round(p1(1,1:2)))+1); %Image indexes start at 1 but CurrentPoint starts at 0,
    p2 = double(uint16(round(p2(1,1:2)))+1); %offset by 1 is to resolve this

    if(p1(1) < p2(1))
        x1 = p1(1);
        x2 = p2(1);
    else
        x2 = p1(1);
        x1 = p2(1);
    end
    if(p1(2) < p2(2))
        y1 = p1(2);
        y2 = p2(2);
    else
        y2 = p1(2);
        y1 = p2(2);
    end

    grab=double(i1(y1:y2,x1:x2));

    bkgn=std(grab(:))/pix_to_pho; %double since adding 2 ROI in localization

    % clear p1 p2 x1 x2 y1 y2 grab infile i1;
end