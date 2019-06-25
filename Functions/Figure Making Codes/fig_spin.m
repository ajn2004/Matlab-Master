% Spin Frames
% A script to rotate and spin a figure by adjusting viewing angles
function fig_spin(fname)
azstrt = 0;
azend = 360;
elstrt = 90;
step = 1;
elend = 0;
% M = [];
clear Mvy

for el = elstrt:-step:elend
    view([azstrt,el]);
    drawnow
    if exist('Mvy')
        Mvy(numel(Mvy)+1) = getframe(gcf);
    else
    Mvy = getframe(gcf);
    end
end

for az = azstrt:step:azend
    view([az,elend])
    drawnow
    if exist('Mvy')
        Mvy(numel(Mvy)+1) = getframe(gcf);
    else
    Mvy = getframe(gcf);
    end
end

for el = elend:step:elstrt
    view([az,el]);
    drawnow
        Mvy(numel(Mvy)+1) = getframe(gcf);

end
movie2gif(Mvy,[fname,'.gif'],'DelayTime',0.02,'LoopCount', Inf);