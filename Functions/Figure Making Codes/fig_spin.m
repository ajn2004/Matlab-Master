% Spin Frames
% A script to rotate and spin a figure by adjusting viewing angles
az = 0;
% el = 32;
% M = [];
clear Mvy

for el = 90:-1:32
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

for el = 32:90
    view([az,el]);
    drawnow
        Mvy(numel(Mvy)+1) = getframe(gcf);

end
movie2gif(Mvy,'neural_gas_hek_cell.gif','DelayTime',0.02,'LoopCount',Inf);