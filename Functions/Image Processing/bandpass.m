% Apply a spatial bandpass filter to an image using gaussians
%
%
function i2 = bandpass(i1,varargin)

if numel(varargin) == 0
    hp = 0.6;
    lp = 1.2;
else
    hp = varargin{1};
    lp = varargin{2};
end

[X,Y] = meshgrid(-5:5);

G1 = exp(-(X.^2+Y.^2)/(2*pi*hp^2));
G2 = abs(1-exp(-(X.^2+Y.^2)/(2*pi*lp^2)));

bp = G1/sum(G1(:)) - G2/sum(G2(:));

i2 = gpu_conv(i1,bp);
i2 = i2.*(i2>0);
end