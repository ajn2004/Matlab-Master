function iprod = bp_subtract(i1, varargin)
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
% se = strel('ball',rball,rball); %structural element, i.e. rolling ball
% i_ball = se.getheight();
% i_hood = se.getnhood();
% kw=10; %kernal width of smoothing function
% [Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);

% xrow = -kw/2:kw/2;
% i_ball= i_ball./sum(sum(i_ball));

% gauss_vec = exp(-2*xrow.*xrow/(rk*rk));
% i_gauss=exp(-2*(Xgs.^2 + Ygs.^2)/(sigma2.^2));
% i_gauss=i_gauss/sum(sum(i_gauss));

% calculate memory requirements
i1 = single(i1);
rows = numel(i1(:,1,1));
cols = numel(i1(1,:,1));
ims = numel(i1(1,1,:));

mem_size = rows*cols*ims*8;

if mem_size < amem / 8             % if memory of the image is small enough, run it through
    [iprod] = bandpass(i1);
else                                                % If memory of image is too small, break it down into chunks and run it through the gpu
    chunks = ceil(mem_size / (amem/10));
    chunkim = ceil(ims / chunks);
    for i = 1:chunks
        if i ~= chunks
            [iprod_temp] = bandpass(i1(:,:,1+(i-1)*chunkim:i*chunkim));
            %             [atemp] = image_neural_2(iprod_temp, Theta1.', Theta2.', numel(iprod_temp(1,1,:)));
            if i ==1
                iprod = iprod_temp;
                %                 a1 = atemp;
                clear iprod_temp
            else
                iprod = cat(3,iprod, iprod_temp);
                %                 a1 = cat(3,a1,atemp);
                clear iprod_temp atemp
            end
        else
            [iprod_temp] = bandpass(i1(:,:,1+(i-1)*chunkim:end));
            %             [atemp] = image_neural_2(iprod_temp, Theta1.', Theta2.', numel(iprod_temp(1,1,:)));
            iprod = cat(3,iprod,iprod_temp);
            %             a1 = cat(3,a1,atemp);
            clear iprod_temp
        end
    end
end
if strcmp(type,'single')
    iprod = single(iprod);
end
disp('Done Subtracting background');