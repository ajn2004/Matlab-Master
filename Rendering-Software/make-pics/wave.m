

[x,y] = meshgrid(-5:0.05:5);
k = 2*pi/1;
R = 4;
r1 = ((x-4).^2 + (y-1).^2).^0.5;
r2 = ((x-4).^2 + (y).^2).^0.5;
r3 = ((x-4).^2 + (y+1).^2).^0.5;
k = 2*pi/1;
w = 2*pi/10;
a = 20;
sig = 50/k;
t = 0;
phi = 0;
while true
t = t + 0.1;
z1 = a*sin(k.*r1-w*t).*exp(-r1.^2/(2*sig^2));
z2 = a*sin(k.*r2-w*t + phi).*exp(-r2.^2/(2*sig^2));
z3 = a*sin(k.*r3-w*t + phi).*exp(-r3.^2/(2*sig^2));
% surf(z)
% surf(z1+z2)
imagesc((z1+z2+z3).^2)
colormap('jet')
axis image
% zlim([-0 a*2])
% view([360/(2*pi)*w/10*t, 30])
drawnow
% M(round(t*10)) = getframe(gcf)
end