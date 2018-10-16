% x-corr example
x = 0:pi/10: 4*pi;
y = sin(x);
% ind = x<0;
% y(ind) = 0;
y1 = sin(x + pi);
% ind = x >0;
% y1(ind) = 0;
plot(x,y,x,y1)
legend('F1',' F2');
[A] = xcorr(y,y1);
for i = 0:2*numel(x)
%     try
%         x1 = x - x(numel(x) - i);
%     catch
        x1 = x + ( i - numel(x)) * pi/10;
%     end
    subplot(1,2,1);
    plot(x,y,x1,y1)
    xlim([-numel(x)*pi/10, 2*numel(x)*pi/10])
    title('visual example')
    subplot(1,2,2)
    plot([-numel(x)*pi/10:pi/10:(i-numel(x))*pi/10]+pi/10,A(1:i+1))
    xlim([-numel(x)*pi/10, numel(x)*pi/10])
    ylim([-20, 15])
    xlabel('Displacement of Functions')
    ylabel('Summed Product of Functions')
    title('Cross Correlation of Functions')
    drawnow
    M(i+1) = getframe(gcf);
end
    