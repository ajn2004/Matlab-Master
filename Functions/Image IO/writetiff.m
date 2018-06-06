function writetiff(i1, fname)
count = 2;
[m,n,p,q] = size(i1);
if ~strcmp(fname(end-3:end),'tif')
    fname = [fname,'.tif'];
end
if q == 1
    imwrite(uint16(i1(:,:,1)),fname)
    if p >1
        while true
            try
                imwrite(uint16(i1(:,:,count)),fname,'WriteMode','append');
                count = count +1;
                if count > numel(i1(1,1,:))
                    break
                end
            catch lsterr
                disp(lsterr);
                lsterr;
            end
        end
    end
else
    imwrite(uint16(i1(:,:,:,1)),fname);
    while true
        try
                imwrite(uint16(i1(:,:,count)),fname,'WriteMode','append');
                count = count +1;
                if count > numel(i1(1,1,:))
                    break
                end
            catch lsterr
                disp(lsterr);
                lsterr;
            end
    end
end