%% Figure for heat explanation
clearvars;

% for k = 0.01:0.01:20
tc = 30;
% tc = k;
hon = 5; % heater on in sec (heater off in sec)
x = 0:0.01:100; % seconds
y1 = x*0; % heater duty cycle
y = exp(-x/tc); % approx response
% y = y1
for i = 0.5:9.5
    ind1 = find(x > i*10 & x <=i*10+5); % heating cycle
    y1(ind1) = 1;
    y(ind1) = (1-exp(-(x(ind1)-i*10)/tc)) + y(ind1(1)-1);
    ind2 = find(x>i*10+5);
    if ~isempty(ind2)
    y(ind2(1):end) = y(ind2(1)-1)*exp(-(x(ind2(1):end)-(i*10+5))/tc);
    end
end
y = y-min(y);
y = y./max(y);
% yg = gausssmooth(y,1,10);
plot(x,y1);
hold on
plot(x,y);
legend('Duty Cycle','Approx Heating','Location','East')
hold off
% title(num2str(k))
ylim([-1,2]);
drawnow
% end

figure
y1 = x*0; % heater duty cycle
% y = exp(-x/tc); % approx response
y = y1;
ind = find(x == 5);
y1(ind:end) = 1;
y(ind:end) = (1-exp(-(x(ind:end) - x(ind))/tc));
plot(x,y1);
figure
plot(x,y1);
hold on
plot(x,y);
