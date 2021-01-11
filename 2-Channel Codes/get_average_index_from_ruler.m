function indexes = get_average_index_from_ruler(ruler, image_stack)
    [m,n,o] = size(image_stack);
    [m,n,r] = size(ruler);
    indexes = zeros(o,1);
    avg_wind = 3; % A window variable so that we average -avg_wind : avg_wind
    %Average image stack
    for i = 1:o
        indexes = i + (-avg_wind:avg_wind);
        indexes(indexes <= 0) = [];
        indexes(indexes > o) = [];
%         disp(indexes)
            average_image_stack(:,:,i) = mean(image_stack(:,:,indexes),3);

    end
    
    for i = 1:o % loop over images
        if round(100*i/o) == 100*i/o
            clc;
            disp(['Indexed ', num2str(100*i/o), '% of file']);
        end
        C = [];
        for j = 1:r
           C(j) = corr2(ruler(:,:,j), average_image_stack(:,:,i));  
        end
        C = gausssmooth(C,5,10);
        try
        indexes(i) = find(C == max(C), 1);
        catch
            indexes(i) = indexes(i-1);
        end

    end
end