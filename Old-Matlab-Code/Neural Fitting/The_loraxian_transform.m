%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Loraxian Transform
%std(
% A transform to take a pixelated image of a psf into the fitted
% localization parameters
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
clear all; clc;
hold on
%% Recreate functional fitting
% [fname, fpath] = uigetfile('*.mat');
% load([fpath,fname]);
load('3d_file.mat');
% for trainper = 0.01:0.01:0.9
count = 1;
clearvars -except x y trainper count 
o = 1;
% x = x.';
for k = 1:o
tic
%% Neural Network Stuff
% Layer information
inputs = numel(x(1,:));
outputs = 1;
hiddens = 37;
it = 1000;

% Fitting Stuff
lambda = 0.1;
trainper = 0.2;
epsint = 1.12;

% data shaping
ind = randperm(numel(x(:,1)));

% indy = find(y == 1.6);
% y(indy) = [];
% x(indy,:) = [];
% xsub = x(ind(1:round(trainper*numel(x(:,1)))),:);
% ysub = y(ind(1:round(trainper*numel(x(:,1)))));
poss = unique(y);
xsub = [];
xf = [];
yfs = [];
ysub =[];
for i = 1:numel(poss)
%     clc;
%     100*i/numel(poss)
    ind = find(y == poss(i));
    xsub = [xsub;(x(ind(1:round(numel(ind)*trainper)),:))];
    ysub = [ysub;(y(ind(1:round(numel(ind)*trainper))))];
    xf = [xf;x(ind(round(numel(ind)*trainper)+1:end),:)];
    yfs = [yfs;y(ind(round(numel(ind)*trainper)+1:end))];
end
% count = 1;
% big test
% for trainper = 0.1:0.01:0.8
% xf = x(ind(round(trainper*numel(x(:,1)))+1:numel(x(:,1))),:);
% yfs = y(ind(round(trainper*numel(x(:,1)))+1:numel(x(:,1))),:);
%         


% for lambda = 1.5:0.1:2.5
% disp([num2str(hiddens),' Hidden Layers'])
% build random thetas
theta1 = rand(hiddens, inputs + 1)*2*epsint - epsint;
theta2 = rand(outputs, hiddens +1)*2*epsint - epsint;
% theta1 = ones(hiddens, inputs + 1);
% theta2 = ones(outputs, hiddens +1);


%feature scaling
[ys,ymi,ysc] = func_feat_scale(ysub);
[theta1, theta2] = func_lorax_neural_teach(theta1, theta2, xsub, ys, lambda, it, 1);


% tic
% trained neural calc
X1 = [ ones(numel(xf(:,1)),1), xf];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a3 = sigmoid(z3);

% disp([num2str(times(count)),'s for computing ', num2str(numel(yfs(:,1))),' molecules'])

X1 = [ones(numel(xsub(:,1)),1), xsub];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a4 = sigmoid(z3);
clear yf yfo yos
for j = 1:1
    % unscale
    yos(:,j) = a3(:,j) * ysc(1,j) + ymi(1,j);
    
    yfo(:,j) = a4(:,j)*ysc(1,j) + ymi(1,j);
    
%     subplot(2,3,j);
%     histogram(yf(:,j)-yfs(:,j));
    
end
% yos = a3;
% yfo = a4;
results{count} =(yos-yfs);
trues{count} = yfs;
tres{count} = (yfo - ysub);
hidden(count) = hiddens;
lambd(count) = lambda;
times(count) = toc;
ajn_wait(times, k, o);
count = count +1;
% end
end

res = [];
tre = [];
for i = 1:count -1
    res= [res; (results{i})];
    tre =[tre; (tres{i})];
end


h = histogram(res,'Normalization','probability');
xlabel('Difference between Truth and Fit in nm');
ylabel('Frequency')
title('Histogram of Results')
drawnow
vals = h.Values;

edg = h.BinEdges;
for i = 1:numel(edg)-1
    cen(i) = 0.5*(edg(i)+edg(i+1));
end
% end
% plot((res))
% hold on
% plot(tre)
% % plot(mean(results) - mean(tres))
% hold off
% histogram(results{:})
% subplot(2,3,1);
% title('Histogram of Difference in xf and truth')
% xlabel('pixels')
% subplot(2,3,2);
% title('Histogram of Difference in yf and truth')
% xlabel('pixels')
% subplot(2,3,3);
% title('Histogram of Difference in n and truth')
% xlabel('photons')
% subplot(2,3,4);
% title('Histogram of Difference in sigx and truth')
% xlabel('pixels')
% subplot(2,3,5);
% title('Histogram of Difference in sigy and truth')
% xlabel('pixels')
% subplot(2,3,6);
% title('Difference in b and truth')
% xlabel('photons')
% save('fitting_thetas.mat','theta1','theta2','ysc','ymi');
% std(results(:,1))*133
% std(results(:,2))*133
% std(results(:,3))
% std(results(:,4))*133
% std(results(:,5))*133
% std(results(:,6))