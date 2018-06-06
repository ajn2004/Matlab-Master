function [iprod] = func_background_sub(i1, rk, rball)
% A function to streamline the gpu background subtraction
% AJN 4/5/17

%% Background Subtraction
rball=5; %radius of rolling ball
amem = 6.0264*10^09;
se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball
i_ball = se.getheight();
i_hood = se.getnhood();
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
% FWHM = 1;
% rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
% xrow = -kw/2:kw/2;
i_ball= i_ball./sum(sum(i_ball));

% gauss_vec = exp(-2*xrow.*xrow/(rk*rk));
i_gauss=exp(-2*kd.*kd./(rk*rk));
i_gauss=i_gauss/sum(sum(i_gauss));

% calculate memory requirements
rows = numel(i1(:,1,1));
cols = numel(i1(1,:,1));
ims = numel(i1(1,1,:));
mem_size = rows*cols*ims*8;
if mem_size < amem / 8             % if memory of the image is small enough, run it through
    [iprod] = image_process(i1,i_gauss, i_ball);
else                                                % If memory of image is too small, break it down into chunks and run it through the gpu
    chunks = ceil(mem_size / (amem/10));
    chunkim = ceil(ims / chunks);
    for i = 1:chunks
        if i ~= chunks
            [iprod_temp] = image_process(i1(:,:,1+(i-1)*chunkim:i*chunkim),i_gauss, i_ball);
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
            [iprod_temp] = image_process(i1(:,:,1+(i-1)*chunkim:end),i_gauss, i_ball);
            %             [atemp] = image_neural_2(iprod_temp, Theta1.', Theta2.', numel(iprod_temp(1,1,:)));
            iprod = cat(3,iprod,iprod_temp);
            %             a1 = cat(3,a1,atemp);
            clear iprod_temp
        end
    end
end
disp('Done Subtracting background');
end