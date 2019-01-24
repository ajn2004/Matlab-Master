% Rolling Ball Efficiency Trials
% close all
% Image Setup
[x,y] = meshgrid(-50:50);
sig1 = 1.5;
sig2 = 3;
sigma2 = 3.9;
rball = 4;
exp1 = exp(-(x.^2/(2*sig1^2)) - y.^2/(4*sig1^2));
exp2 = exp(-(x.^2+y.^2)/(2*sig2^2));

im = exp1+exp2;
% rim = rollingball(im,1,4);
% rball =
se = strel('ball',rball,rball); %structural element, i.e. rolling ball
i_ball = single(se.getheight());
i_hood = single(se.getnhood());
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
i_gauss=exp(-2*(Xgs.^2 + Ygs.^2)/(sigma2.^2));
i_gauss=i_gauss/sum(sum(i_gauss));

imc = conv2(im,i_gauss,'Same');
imo = imopen(imc,se);
rim = im - imo;

subplot(3,2,1);
surf(exp1);
title('Gauss 1');
% axis image
subplot(3,2,5);
surf(exp2);
title('Gauss 2');
% axis image
subplot(3,2,3);
surf(im);
% axis image
title('Raw Image')
subplot(3,2,2);
surf(rim);
view([-11 0])
zlim([-1,2])
title('Background Subtracted')
subplot(3,2,4);
surf(rim-exp1);
zlim([-1,2])
view([-11 0])
title('BS - exp1')
subplot(3,2,6);
surf(rim-exp2);
zlim([-1,2])
title('BS-exp2');
view([-11 0])

% axis image


