function [iprod, ibkgn] = rollingball(i1, varargin)
%ROLLINGBALL a function to perform gpu rolling ball subtraction on an image
%   A = ROLLINGBALL(IM1) will return a matrix A the same size as IM1 that
%   has been background subtracted using the rolling ball method with
%   defaults of rball = 5 and sigma = 2
%   
%   A = ROLLINGBALL(IM1,RBALL,SIGMA) will return a matrix A the same size as IM1 that
%   has been background subtracted using the rolling ball method with the
%   specified variables
type = class(i1);
numvar = length(varargin);
if numvar > 2
    error('Too many variables');
elseif numvar <2
    rball = 5;
    sigma2 = 4;
elseif numvar == 2
    rball = varargin{1};
    sigma2 = varargin{2};
end
gpud = gpuDevice;
amem = gpud.AvailableMemory;
se = strel('ball',rball,rball); %structural element, i.e. rolling ball
i_ball = single(se.getheight());
i_hood = single(se.getnhood());
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);


i_ball= i_ball./sum(sum(i_ball));


i_gauss=exp(-2*(Xgs.^2 + Ygs.^2)/(sigma2.^2));
i_gauss=single(i_gauss/sum(sum(i_gauss)));

% calculate memory requirements
i1 = single(i1);
[rows, cols, ims] = size(i1);
mem_size = rows*cols*ims*8;

if mem_size < amem / 8             % if memory of the image is small enough, run it through
    [iprod,ibkgn] = image_process(i1,i_gauss, i_ball);
else                                                % If memory of image is too small, break it down into chunks and run it through the gpu
    chunks = ceil(mem_size / (amem/10));
    chunkim = ceil(ims / chunks);
    for i = 1:chunks
        if i ~= chunks
            [iprod_temp,bkgn_temp] = image_process(i1(:,:,1+(i-1)*chunkim:i*chunkim),i_gauss, i_ball);
            
            if i ==1
                iprod = iprod_temp;
                ibkgn = bkgn_temp;
            
                clear iprod_temp
            else
                iprod = cat(3,iprod, iprod_temp);
                ibkgn = cat(3,ibkgn, ibkgn_temp);
            
                clear iprod_temp atemp
            end
        else
            [iprod_temp,bkgn_temp] = image_process(i1(:,:,1+(i-1)*chunkim:end),i_gauss, i_ball);
            
            iprod = cat(3,iprod,iprod_temp);
            ibkgn = cat(3,ibkgn, ibkgn_temp);
            
            clear iprod_temp
        end
    end
end
if strcmp(type,'single')
    iprod = single(iprod);
    ibkgn = single(ibkgn);
end
disp('Done Subtracting background');