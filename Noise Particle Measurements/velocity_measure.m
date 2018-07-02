clearvars; close all; clc;

[path, bright] = dlinetrace();

vel = [0];
try
while true
[x,y] = ginput(2);

dist = 0.133*sum((path(round(y(1)),:) - path(round(y(2)),:)).^2)^0.5;
vel(numel(vel)+1) = dist/(0.0409*abs(round(x(1))-round(x(2))));
end
catch lsterr
    vel(1) = [];
end
disp(['Calculated Velocity is ', num2str(mean(vel)),'+/-',num2str(std(vel)/numel(vel)^0.5),'um/s over ', num2str(numel(vel)), ' measurements']);