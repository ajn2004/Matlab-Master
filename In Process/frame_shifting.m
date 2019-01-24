ts = 50:25:949; % list the frames for the stimuli
% trying to find the frames of photobleach by looking at the difference i
% frames after the stimulus and taking the average value and finding a
% regular drop
clear xfluor
wind = -4:20;
for i = 1:numel(ts)
    xfluor(:,i) = mfluor(ts(i)+wind)-mean(mfluor(ts(i)-4:ts(i)));
end
% errorbar(1:i,mi,ste)
plot(xfluor,'Color',[0.8,0.8,0.8])
hold on
plot(mean(xfluor,2),'Color',[0 0 0])
hold off
