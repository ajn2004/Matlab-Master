function i2 = gpu_rball(i1,rball,g_sig)
if nargin<3
    rball = 6;
    g_sig = 5;
end
[m,n,o] = size(i1);
pixels = m*n*o; % number of pixels in image
i1 = double(i1); % ensure the image is of type double
maxpix = 224*228*1000; % empirical memory maximum before violating gpu
i2 = [];

% Processing elements
se = strel('ball',rball,rball); %structural element, i.e. rolling ball
i_ball = double(se.getheight());
i_hood = double(se.getnhood());
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
i_gauss=exp(-2*(Xgs.^2 + Ygs.^2)/(g_sig.^2));
i_gauss=i_gauss/sum(sum(i_gauss));

if pixels > maxpix % if there are more pixels than the memory limit
    chunks = ceil(pixels/maxpix);  % break the process up into a number of chunks
    imchunks = round(o/chunks);
    for i = 1:chunks + 1
        im2 = [];
        try
           [im2, imo] = gpu_rolling_ball_11(i1(:,:,imchunks*(i-1)+1:imchunks*i),i_gauss,i_ball);
        catch
            if imchunks*(i-1) +1 < o
                [im2, imo] = gpu_rolling_ball_11(i1(:,:,imchunks*(i-1)+1:end),i_gauss,i_ball);
            end
        end
        i2 = cat(3,i2,im2);
    end
    
else
    [i2,imo] = gpu_rolling_ball_11(i1, i_gauss, i_ball);
end
