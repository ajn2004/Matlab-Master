% Apply a gaussian blur
%
%
function i2 = imblur(i1,varargin)

if numel(varargin) == 0
    sig = 3;
elseif nargin == 1
    sig = varargin{1};
else
    sig = varargin{1};
end

[X,Y] = meshgrid(-5:5);

G1 = exp(-(X.^2+Y.^2)/(2*sig^2));

bp = G1/sum(G1(:));
i2 = image_conv(i1,bp);
% i2 = i2.*(i2>0);
end