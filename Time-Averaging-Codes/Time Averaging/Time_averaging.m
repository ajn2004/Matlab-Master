% Trial Averaging
%
% Averages trials together to get a final graph
%

[fname, fpath] = uigetfile('*.mat');
stim = 50;
cd(fpath);
finfo = dir('*.mat');
singles = [];
doubles = [];
for i = 1:numel(finfo)
    load(finfo(i).name);
    nstr = finfo(i).name;
    singles = [singles;fluor];
end

figure
ave_sing = mean(singles,1);
x = 1:numel(ave_sing);
gas = gausssmooth(ave_sing(:), 1.4,10);
plot(x*0.01,gas);
% legend('Signal');
% legend('Double','Single')
title('Averaged response to stimulus')
ylabel('dF / F')
xlabel('Time in (s)');
asings = ave_sing(:);
hold on
for i = 1:5
    plot(0.01*[(i-1)*10 + 50, (i-1)*10 + 50], [-0.001, 0.006],'r');
end

plot(0.01*[40, 55], [mean(gas(44:53)), mean(gas(44:53))],'k','LineWidth',2);
plot(0.01*[55, 70], [mean(gas(57:70)), mean(gas(57:70))],'k','LineWidth',2);
plot(0.01*[70, 80], [mean(gas(74:78)), mean(gas(74:78))],'k','LineWidth',2);
plot(0.01*[75, 90], [mean(gas(79:88)), mean(gas(79:88))],'k','LineWidth',2);
plot(0.01*[85, 100], [mean(gas(90:98)), mean(gas(90:98))],'k','LineWidth',2);
hold off

% adubs = ave_dub(:);

% sing_fit = asings(stim+2:end);
% dub_fit = adubs(stim+2:end);

% % sfit = polyfit([stim+2:numel(ave_dub)].',sing_fit,1);
% % dfit = polyfit([stim+2:numel(ave_dub)].',dub_fit,1);
% x = 1:numel(asings);
% sy = sfit(1)*x + sfit(2);
% % dy = dfit(1)*x + dfit(2);
% plot(100*sy);
% % plot(100*dy);
% hold off
% % separ = dfit(2)/sfit(2);
% disp(['The ratio of y-intercepts in the lines fitted to the linear sections is ' , num2str(separ)]);

% rats = dub_fit./sing_fit;
% mean(rats)