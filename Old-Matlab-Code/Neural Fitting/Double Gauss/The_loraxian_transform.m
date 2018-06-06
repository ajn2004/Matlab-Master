%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Loraxian Transform
%std(
% A transform to take a pixelated image of a psf into the fitted
% localization parameters
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Recreate functional fitting
[fname, fpath] = uigetfile('*.mat');
load([fpath,fname]);


%% Neural Network Stuff
% Layer information
inputs = 81;
outputs = 1;
hiddens = 25;
it = 5000;

% Fitting Stuff
lambda = 0;
trainper = 0.4;
epsint = 0.12;

% data shaping
ind = randperm(numel(i2(:,1)));

xsub = i2(ind(1:round(trainper*numel(i2(:,1)))),:);
ysub = truth(ind(1:round(trainper*numel(i2(:,1)))),:);

% big test
xf = i2(round(trainper*numel(i2(:,1)))+1:numel(i2(:,1)),:);
yfs = truth(round(trainper*numel(i2(:,1)))+1:numel(i2(:,1)),:);
        
count = 1;
% for hiddens = [15, 17, 19, 21, 23, 25, 27]
disp([num2str(hiddens),' Hidden Layers'])
% build random thetas
theta1 = rand(hiddens, inputs + 1)*2*epsint - epsint;
theta2 = rand(outputs, hiddens +1)*2*epsint - epsint;



%feature scaling
[ys,ymi,ysc] = func_feat_scale(ysub);
[theta1, theta2] = func_lorax_neural_teach(theta1, theta2, xsub, ys, lambda, it, 1);


tic
% trained neural calc
X1 = [ ones(numel(xf(:,1)),1), xf];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a3 = sigmoid(z3);
times(count) = toc;
disp([num2str(times(count)),'s for computing ', num2str(numel(yfs(:,1))),' molecules'])

for j = 1:numel(a3(1,:))
    % unscale
    yf(:,j) = a3(:,j) * ysc(1,j) + ymi(1,j);
    
    
    
%     subplot(2,3,j);histogram(yf(:,j)-yfs(:,j));
    
end

results(:,:,count) =(yf-yfs);
hidden(count) = hiddens;
count = count +1;
tend = mean(round(results) == 0 ,1)
plot(hidden, tend(:));
drawnow;
% end
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