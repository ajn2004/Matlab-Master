clearvars
close all
pixw = 20;
x = zeros(400,1);
x(100:199) = 1;
x(200:299) = -1;
x(300:400) = 3;

% kernel construction
xx = -pixw:pixw; % create x axis
sx = pi*xx/pixw; % convert to -pi - pi to ensure area under sine is 0
sig = pixw;     % assign sigma value
gs = -sin(sx).*gaussian(sig,pixw); % calculate gaussian multiplied by sine curve for edge filter
edges = conv(x,gs,'same'); % convolve F with gauss-sine kernel and return the result

for i = 1:400
    plot(x)
    hold on
    plot(i+xx,-10*gs,'k')
    plot(1:i,edges(1:i),'r')
    hold off
    xlim([1,400])
    title('Example Edge Detection');
    xlabel('Arbitrary X')
    ylabel('Arbitrary Y')
    drawnow
    M(i) = getframe(gcf);
end
    


