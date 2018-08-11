function i2 = xdrift(i1)
amem = 147*146*60; %safe amount of safe pixels for driftgpu calc
[m,n,o] = size(i1);
nmem = m*n*o;

if nmem <= amem
    i2 = driftgpu(i1,1);
else
    i2 = [];
    chunks = round((nmem -1)/amem +1);
    imchk = ceil(o/chunks);
    for i = 1:chunks
        try
            im2 = driftgpu(i1(:,:,(i-1)*imchk+1:i*imchk),10);
            i2 = cat(3,i2,im2);
        catch
            im2 = driftgpu(i1(:,:,i*imchk:end),10);
            i2 = cat(3,i2,im2);
        end
    end
    
end