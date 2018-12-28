% TLR Script
close all;
clearvars
s = (0.61*0.51/1.3)/2;
a = 0.08;
% N = 50:50000;
N = 200;
snr = 20:0.1:100;
alpha = 20;
% i=100;
% for i = 0:5 % loop over number of vesicles in a buton
    b = 0.1*N; % background is equal to the dark rate of phluorin
    N  = 1500;
    c = jet(numel(snr));
%     for i  = 1:numel(snr)
    lp2 = (s^2+a^2/12)./N + (8*pi*s^4)./(a^2*snr.^2);
    lp = lp2.^0.5;
%     plot(N,lp*1000,'Color',c(i,:))
plot(snr,lp*1000)
%     plot(snr,lp*1000,'Color',c(i,:))
    hold on
%     end
% end
plot([snr(1),snr(end)],[0.015,0.015]*1000,'r')
% plot([N(1),N(end)],[0.02,0.02]*1000,'r')
xlim([min(snr),max(snr)])
% xlim([min(N),max(N)])
ylim([min(lp)*1000,1000*max(lp)]);
ylim([9, 45]);
% legend('No Background','1 Molecule','2 Molecules','3 Molecules','4 Molecules','5 Molecules')
set(gca,'YScale','log');
% set(gca,'XScale','log');
xlabel('SNR');
ylabel('Loc Precision nm')
grid on
hold off