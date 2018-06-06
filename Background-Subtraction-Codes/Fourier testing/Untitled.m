% timing tester

for i = 2:7
    i1 = zeros(10^i,100*100);
    tic
    fi = fft(i1,10^i,1);
    fi(1:20,:) = 0;
    fi(end-21:end,:) = 0;
    i2 = ifft(fi,10^i,1);
    clear fi
    tc(i-1) = toc
    
    tic
    iprod = cgpufourier(i1,20);
    tg(i-1) = toc
end
x = 2:7;
plot(10.^x,tc);
hold on
plot(10.^x,tg,'r');
legend('CPU','GPU');
hold off
    