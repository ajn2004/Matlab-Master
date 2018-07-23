% focus lock sim
x = 45:0.1:55;
y = 45:0.1:55;
close all
dx = 0.05;
dy = 0.05;

for j = 1:numel(x)
    for k = 1:numel(y)
        nfree(j,k) = func_lock(x(j),y(k));
        nfre(j,k) = x(j)/y(k);
        clear nse nsse
        for i = 1:10000
            nx = x(j)*randn*dx;
            ny = y(j)*randn*dy;
            
            nse(i) = nfree(j,k) - func_lock(x(j) + nx,y(k) + ny);
            nsse(i) = nfre(j,k) - (x(j) + nx)/(y(k) + ny);
        end
        std_noiser(j,k) = std(nsse);
        std_noise(j,k) = std(nse);
    end
    j
end
surf(std_noise);
title('Difference over sum')
figure
surf(std_noiser);
title('Straight Ratio')

figure
surf(std_noise./std_noiser)
title('Ratio of methods')