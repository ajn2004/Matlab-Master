%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic Curve Fitting
%
% A transform to take a pixelated image of a psf into the localized
% coordinates
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; clc; close all;
% hold off
%% Use this section to create your data set 
% load('composit.mat','x','drifts')
clear yf theta1 theta2 xsub ysub

%% 
% load('3d_data.mat');
% ind = absz < 1000;
x=1:0.01:10;
x = x(:);
y = [x+2,mod(x,3)];
xsub = x;
ysub = y;
for i = 1:numel(ysub(1,:))
plot(x,y(:,i),'.')
hold on
end

%% Neural Network Stuff 
% Layer information
inputs = numel(xsub(1,:));
outputs = numel(ysub(1,:));
hiddens = 30;
it = 5000;

% Fitting Stuff
lambda = 0.0;
trainper = 1;
epsint = 0.12;

% build random thetas
theta1 = rand(hiddens, inputs + 1)*2*epsint - epsint;
theta2 = rand(outputs, hiddens +1)*2*epsint - epsint;

% Feature scaling
for i = 1:numel(ysub(1,:))
    ymi(i) = min(ysub(:,i));
    ysc(i) = (max(ysub(:,i)) - min(ysub(:,i)));
    ys(:,i) = (ysub(:,i) - ymi(i))/ysc(i);
end

[theta1, theta2] = func_lorax_neural_teach(theta1, theta2, xsub, ys, lambda, it, trainper);

% trained neural calc
% x = (1:100).';
a3 = func_eval_NN(x,theta1,theta2);
% x = x(:);
% X1 = [ ones(numel(x(:,1)),1), x];
% z2 = X1*theta1.';
% a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
% z3 = a2*theta2.';
% a3 = sigmoid(z3);

% unscale
for i = 1:numel(ysub(1,:))
yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
end

for i = 1:numel(ysub(1,:))
hold on
plot(x,yf(:,i),'.');
end
legend('Data',' Neural Network');