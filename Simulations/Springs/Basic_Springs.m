% Spring simulations
clearvars; close all; clc;

x1 = -1;
x2 = 1;
k1 = 1;
k2 = 1;
p = 0.5;
El = 300;
E = 0;
dt = 0.001;
count = 1;
while abs(E) > -1
    El = E;
    E1 = k1*(x1-p)^2;
    E2 = k2*(x2-p)^2;
    E = E1 + E2;
    F1 = -k1*(x1-p);
    F2 = -k2*(x2-p);
    F = F1 + F2;
    dx = dt*2*F;
    p = p - dx;
    plot(0,0,'rx','MarkerSize',10)
    hold on
    plot([x1, p, x2], [0 0 0], '.')
    axis equal
    xlabel(num2str(E))
    ylim([-1 1])
    xlim([x1, x2]*1.1)
    drawnow
    hold off
    M(count) = getframe(gcf);
    count = count +1;
end
    