% putstims
close all
% glfo = gausssmooth(mfluor,20,10);
glfo = mfluor;
% h = [-1 -1 -1 - 1 0 1 1 1 1];
pixw = 5;
wind = -pixw:pixw;
sind = (wind/pixw)*pi;
h = sin(sind).*gaussian(pixw,pixw);
h1 = sigmoid(wind)-0.5;
cflo = conv(glfo,-h,'same');
cflo1 = conv(glfo,-h1,'same');
mxfl = [min(cflo1),max(cflo1)];
% mxfl = [min(mfluor),max(mfluor)];
stim = 48;
% plot(cflo);
% hold on
% plot(gausssmooth(cflo1,5,10))
% hold on
plot(cflo1)

% plot(mfluor)
% gflur = gausssmooth(mfluor,500,200);
% nfluor = mfluor - gflur;
% plot(nfluor)
% mxfl = [min(nfluor),max(nfluor)];
hold on
for i = 1:floor(numel(mfluor)/stim)
    plot([stim,stim]*i,mxfl,'r');
    stims(i) = stim*(i);
end
hold off
% ylim([min(cflo(stim:end-stim)) max(cflo(stim:end-stim))])
figure
nel = numel(stims(2:end-1));
histogram(cflo(stims(2:end-1)),2*round(nel^(2/3)))
hold on
histogram(cflo1(stims(2:end-1)),2*round(nel^(2/3)))

figure
plot(mfluor);
mxfl = [min(mfluor),max(mfluor)];
hold on
for i = 1:floor(numel(mfluor)/stim)
    plot([stim,stim]*i,mxfl,'r');
    stims(i) = stim*(i);
end
hold off