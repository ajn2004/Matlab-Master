function indexes = get_index_from_ruler(ruler, image_stack)
    [m,n,o] = size(image_stack);
    [m,n,r] = size(ruler);
    indexes = zeros(o,1);
    for i = 1:o % loop over images
        if round(100*i/o) == 100*i/o
            clc;
            disp(['Indexed ', num2str(100*i/o), '% of file']);
        end
        C = [];
        for j = 1:r
           C(j) = corr2(ruler(:,:,j), image_stack(:,:,i));  
        end
        C = gausssmooth(C,5,10);
        indexes(i) = find(C == max(C));

    end
end