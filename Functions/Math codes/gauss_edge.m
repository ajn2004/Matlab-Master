function edges = gauss_edge(F,pixw)
% this function will take in an array F and a pixel window pixw to perform
% and edge filtration using a gaussian multiplied by a sin wave

% kernel construction
x = -pixw:pixw; % create x axis
sx = pi*x/pixw; % convert to -pi - pi to ensure area under sine is 0
sig = pixw;     % assign sigma value
gs = -sin(sx).*gaussian(sig,pixw); % calculate gaussian multiplied by sine curve for edge filter
edges = conv(F,gs,'same'); % convolve F with gauss-sine kernel and return the result