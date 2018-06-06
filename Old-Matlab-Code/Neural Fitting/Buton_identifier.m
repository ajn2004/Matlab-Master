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
% load([fpath,fname]);
cd(fpath);
%% Data loading and reshaping
i2 = [];
truth = [];

finfo = dir('*.mat');
for i = 1:numel(finfo)
    load(finfo(i).name);
    i2 = [i2;xf];
    truth = [truth;yf(:,1)];
    clear xf yf
end

%% Neural Network Stuff
% Layer information
inputs = 21;
outputs = 1;
hiddens = 10;
it = 10000;

% Fitting Stuff
lambda = 2;
trainper = 0.6;
epsint = 0.12;

% data shaping
ind = randperm(numel(i2(:,1)));

xsub = i2(ind(1:round(trainper*numel(i2(:,1)))),:);
ysub = truth(ind(1:round(trainper*numel(i2(:,1)))),:);

% big test
xf = i2(round(trainper*numel(i2(:,1)))+1:numel(i2(:,1)),:);
yfs = truth(round(trainper*numel(i2(:,1)))+1:numel(i2(:,1)),:);
        
count = 1;
% for hiddens = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200,225,250]
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

for j = 1
    % unscale
    yf(:,j) = a3(:,j) * ysc(1,j) + ymi(1,j);
    
    
    
    subplot(2,3,j);histogram(yf(:,j)-yfs(:,j));
    
end

results(:,:,count) =(round(yf)==yfs);
hidden(count) = hiddens;
% count = count +1
% end
mean(results)
std(results(:,1))
save('buton_thetas.mat','theta1','theta2');
% std(results(:,2))*133
% std(results(:,3))
% std(results(:,4))*133
% std(results(:,5))*133
% std(results(:,6))