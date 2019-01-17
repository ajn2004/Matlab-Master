% Apply a spatial bandpass filter to an image using gaussians
%
%
function i2 = bandpass(i1,varargin)

if numel(varargin) == 0
    hp = 1.5;
    lp = 2.5;
else
    hp = varargin{1};
    lp = varargin{2};
end

[X,Y] = meshgrid(-13:13);

G1 = exp(-(X.^2+Y.^2)/(2*hp^2));
G2 = exp(-(X.^2+Y.^2)/(2*lp^2));

% bp = G1/sum(G1(:)) - G2/sum(G2(:));

bp = -1/(pi*lp^4)*(1-(X.^2 + Y.^2)/(2*lp^2)).*G2;
bp = -1/(pi*hp^4)*(1-(X.^2 + Y.^2)/(2*hp^2)).*G1;

% i2 = gpu_conv(i1,bp);
[m,n,o] = size(i1);
for i = 1:o
    i2(:,:,i) = conv2(i1(:,:,i),-bp,'same');
end
i2 = i2.*(i2>0);
end