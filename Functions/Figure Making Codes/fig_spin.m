% Spin Frames
% A script to rotate and spin a figure by adjusting viewing angles
function fig_spin(fname)
azstrt = 20;
azend = 360;
elstrt = 0;
step = 1;
elend = 0;
% M = [];
clear Mvy

for el = elstrt:-step:elend
    view([azstrt,0]);
    drawnow
    if exist('Mvy')
        Mvy(numel(Mvy)+1) = getframe(gcf);
    else
    Mvy = getframe(gcf);
    end
end

for az = azstrt:step:azend
    view([az,0])
    drawnow
    if exist('Mvy')
        Mvy(numel(Mvy)+1) = getframe(gcf);
    else
    Mvy = getframe(gcf);
    end
end

for el = elend:step:elstrt
    view([az,0]);
    drawnow
        Mvy(numel(Mvy)+1) = getframe(gcf);

end
movie2gif(Mvy,[fname,'.gif'],'DelayTime',0.02,'LoopCount', Inf);