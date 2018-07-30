%% Figure for heat explanation

x = 1:100;
x = x*0;
for i = 1:10
    x(i*10:i*10+4) = 1;
end

plot(x)
ylim([-1,2]);
