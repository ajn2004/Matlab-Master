% Spin Frames
% A script to rotate and spin a figure by adjusting viewing angles
az = 0;
% el = 90;
% M = [];
clear Mvy

for el = 90:-1:20
    view([az,el]);
    drawnow
    if exist('Mvy')
        Mvy(numel(Mvy)+1) = getframe(gcf);
    else
    Mvy = getframe(gcf);
    end
end

for az = 0:360
    view([az,el])
    drawnow
    Mvy(numel(Mvy)+1) = getframe(gcf);
end

for el = 20:90
    view([az,el]);
    drawnow
        Mvy(numel(Mvy)+1) = getframe(gcf);

end
movie2gif(Mvy,'3um_hek_cell.gif','DelayTime',0.02,'LoopCount',Inf);