%% Localization with stimuli
% a code to show which localizations correspond to stimuli and where they
% go afterwards. It requires understanding frame information. If it was
% created from a cut_up_data set additional.mat should have the true frames
% in the form of fms. This will probably serve as the basis for future code
clearvars; close all; clc;
[dname, dpath] = uigetfile('*_tol*');
[aname, apath] = uigetfile('*additional*');
load([dpath, dname]);
load([apath, aname]);
Points_diag;
close all
cfms = fms(x); % corrected frames
mfms = mod((cfms - 100),64).'; % modded frames, subtract the prescans then mod 64 which is the stim rate
xf = xf - mean(xf);
yf = yf - mean(yf);
zf = zf - mean(zf);
plot3(xf,yf,zf,'.b');
ind = mfms <= 1;
hold on
plot3(xf(ind),yf(ind),zf(ind),'.r','MarkerSize',20);
%% Plotting 'Trajectories'
used = [];
for i = 1:numel(zf)
    if ~ismember(i,used)
        [ind] = find(cfms(i) <= cfms + 10 & cfms(i) >=cfms);
        used = [used;ind.'];
        scatter3(xf(ind),yf(ind),zf(ind),[],get_color(numel(ind)));
        plot3(xf(ind),yf(ind),zf(ind))
    end
end


hold off
legend('All locs','Early Locs');
axis equal

title('Rotated Bouton');
xlabel('Position [um]');
ylabel('Position [um]');
zlabel('Position [um]');
view([0,90]);
M(1) = getframe(gcf);
eles=(90:-1:0);

for i = 1:numel(eles)
    view([0,eles(i)]);
    M(numel(M)+1) = getframe(gcf);
end

azs = 1:360;
for i = 1:numel(azs)
    view([azs(i),0]);
    M(numel(M)+1) = getframe(gcf);
end