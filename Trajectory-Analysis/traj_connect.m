function [traj, foll] = traj_connect(foll, traj, strt)
    if strt ~= 0
        traj = [traj;strt];
        [traj, foll] = traj_connect(foll, traj, foll(strt));
        foll(foll == strt) = 0;
    end
end