function i2 = padthis(i1,pad)
% will add an apron of zeros to an image i1 the size of pad
i2 = [];
[m,n] = size(i1);
vertp = zeros(pad,n);
horzp = zeros(2*pad + m,pad);
i2 = [horzp,[vertp;i1;vertp],horzp];
end