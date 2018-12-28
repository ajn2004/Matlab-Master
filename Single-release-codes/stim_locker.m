%Stim_lock
fps = 25;
fstim = 50;
close all
stims = fstim:fps:numel(mfluor)-fps;
histogram(mfluor(stims+1) - mfluor(stims),20);
thsh = 900
ind = mfluor(stims+1) - mfluor(stims) > thsh;
wind = -5:fps-5;
msub = wind*0;
nsub = msub;
for i = 1:numel(stims)
    fsub = mfluor(stims(i)+wind);
    % normalized
    mb = mean(fsub(1:5));
    fsub = fsub - mb;
    fsub = fsub/max(fsub);
%     fsub = fsub/mb;
    if ind(i) == 1
        msub = msub + fsub.';
        psubs(i,:) = fsub;
    else
        nsub = nsub + fsub;
        nsubs(i,:) = fsub;
    end
end

msub = msub / sum(ind);
nsub = nsub / sum(1-ind);
figure
plot(wind*tex,psubs(ind,:),'Color',[0.8,0.8,0.8])
hold on
plot((wind)*tex,msub,'b')
ylim([-1, 1.1])
xlabel('Time in Seconds');
ylabel('Normalized Fluorescence')
title('Time Locked Single Stimulations')
hold off
figure
plot(wind*tex,nsubs(logical(1-ind),:),'Color',[0.8,0.8,0.8])
hold on
plot((wind)*tex,nsub,'r')
hold off
xlabel('Time in Seconds');
ylabel('Normalized Fluorescence')
title('Time Locked Failures')