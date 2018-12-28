function i2 = roball(i1,sigma2,rball)
i2 = i1*0;
[m,n,o] = size(i2);
se = strel('ball',rball,rball); %structural element, i.e. rolling ball
i_ball = single(se.getheight());
i_hood = single(se.getnhood());
kw=10; %kernal width of smoothing function
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
i_gauss=exp(-2*(Xgs.^2 + Ygs.^2)/(sigma2.^2));
i_gauss=i_gauss/sum(sum(i_gauss));

for i = 1:o
    imc = conv2(i1(:,:,i),i_gauss,'Same');
    imo = imopen(imc,se);
    i2(:,:,i) = i1(:,:,i)-imo;
end

i2 = i2.*(i2>0);