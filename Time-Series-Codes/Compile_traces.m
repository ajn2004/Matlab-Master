clearvars
close all
files = dir('*mat');
hz = 45.57;
afluor = [];
sfluor = [];

%% Calculating gaussian window for smoothing
ft = 1/hz;
stim_rate = 1;
stim_p = 1/stim_rate;
fps = stim_p/ft;
gs = 30*fps;
for i = 1:numel(files)
    load(files(i).name);
    if ~isempty(afluor)
        cutf = min([numel(afluor(:,1)),numel(ifluor(:,1))]);
        afluor = [afluor(1:cutf,:),ifluor(1:cutf,:)];
        for j = 1:numel(ifluor(1,:))
            sfluor = [sfluor(1:cutf,:),ifluor(1:cutf,j) - gausssmooth(ifluor(1:cutf,j),gs,round(3*gs))];
        end
    else
        afluor = ifluor;
        for j = 1:numel(ifluor(1,:))
            sfluor = [sfluor,ifluor(:,j) - gausssmooth(ifluor(:,j),gs,round(3*gs))];
        end
    end
    
    clear ifluor
end

mt = (1:numel(sfluor(:,1)))*ft;
mafluor = mean(afluor,2);
msfluor = mean(sfluor,2);
showfft(msfluor,mt);
figure
plot(mt,mafluor)
figure
plot(mt,msfluor);
title('Subtracted')