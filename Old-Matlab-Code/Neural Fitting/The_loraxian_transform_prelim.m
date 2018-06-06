%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic Curve Fitting
%
% A transform to take a pixelated image of a psf into the localized
% coordinates
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Recreate functional fitting
load('composit.mat','x','drifts')

xsub = x;
ysub = drifts;
for i = 1:numel(ysub(1,:))
plot(xsub,ysub(:,i))
hold on
end

%% Neural Network Stuff 
% Layer information
inputs = numel(xsub(1,:));
outputs = numel(ysub(1,:));
hiddens = 20;
it = 1000;

% Fitting Stuff
lambda = 0.00001;
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
x = 1:10000;
x = x(:);
X1 = [ ones(numel(x(:,1)),1), x];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a3 = sigmoid(z3);

% unscale
for i = 1:numel(ysub(1,:))
yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
end

hold on
plot(x,yf,'.r');
legend('Data',' Neural Network');