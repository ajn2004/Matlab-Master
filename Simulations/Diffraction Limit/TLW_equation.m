% TLR Script
clearvars
s = (0.61*0.51/1.3)/2;
a = 0.08;
N = 50:50000;
alpha = 20;
for i = 0:10:300 % loop over number of vesicles in a buton
    b = i*N/alpha; % background is equal to the dark rate of phluorin
    lp2 = (s^2+a^2/12)./N + (8*pi*s^4*b)./(a^2*N.^2);
    lp = lp2.^0.5;
    plot(N,lp*1000)
    hold on
    
end
plot([N(1),N(end)],[0.02,0.02]*1000,'r')
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('Photons Detected');
ylabel('Loc Precision nm')
grid on
hold off