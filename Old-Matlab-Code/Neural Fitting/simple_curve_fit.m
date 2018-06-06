% simple curve fitting
clearvars; clc; close all;

[fname, fpath] = uigetfile('*.mat');
load([fpath,fname]);
tp = 0.6;
% ysub = x.';
% xsub = [sigx.',sigy.'];
ind = randperm(numel(ysub(:,1)));
x = xsub(ind(1:round(tp*numel(xsub(:,1)))),:);
y = ysub(ind(1:round(tp*numel(xsub(:,1)))),:);

[t1, t2,ysc,ymi] = get_them_thetas(x,y,12,1,0.12);
xo = xsub(ind(1+round(tp*numel(xsub(:,1)))):end,:);
a3 = func_eval_NN(xo,t1,t2);
% unscale
yf = [];
for i = 1:numel(ysub(1,:))
    yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
end
subplot(2,1,1);
plot(ysub,xsub(:,1),'.');
hold on
plot(ysub,xsub(:,2),'.');
[yp,  I]  = sort(yf);
plot(yp,xo(I,1));
plot(yp,xo(I,2));
hold off
% legend('sigx','sigy','fitx','fity');
subplot(2,1,2);
histogram(y)
save('n_thets.mat','t1','t2');

