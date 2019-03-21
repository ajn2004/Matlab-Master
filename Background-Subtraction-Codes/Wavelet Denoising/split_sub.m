function im2 = split_sub(i1,i2)
[m,n,o] = size(i1);
pix = m*n*o;
maxpix = 224*228*5000;
im2 = [];
if pix<=maxpix
    im2 = wave_sub(i1,i2,2);
else
    chunks = ceil(pix/maxpix);
    imch = round(o/chunks);
    for i = 1:chunks +1
        imo = [];
        try
            imo = wave_sub(i1(:,:,imch*(i-1)+1:imch*i),i2(:,:,imch*(i-1)+1:imch*i),2);
        catch
            if imch*(i-1) + 1 < o
                imo = wave_sub(i1(:,:,imch*(i-1)+1:end),i2(:,:,imch*(i-1)+1:end),2);
            end
        end
        im2 = cat(3,im2,imo);
    end
            
end