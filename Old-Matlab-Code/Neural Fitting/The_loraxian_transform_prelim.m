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
% load('composit.mat','x','drifts')
x = -2*pi:0.01:2*pi;
% y = abs(sin(x)./(x));%+0.15*randn(1,numel(x));
y = sin(x)
% y = x.^2 + randn(1,numel(x))*2;

inds = randperm(numel(x));
xsub = x(inds(1:round(0.15*numel(x)))).';
ysub = y(inds(1:round(0.15*numel(x)))).';
% xsub = x;
% ysub = drifts;
for i = 1:numel(ysub(1,:))
plot(xsub,ysub(:,i),'.')
hold on
end
hold off
%% Neural Network Stuff 
% Layer information
inputs = numel(xsub(1,:));
outputs = numel(ysub(1,:));
hiddens = 15;
it = 5000;

% Fitting Stuff
lambda = 00.01;
trainper = 1;
epsint = 0.12;
its = [1,5,7,9, 10,15,20,25,30,35,40,50,100,500,1000];
for k = 1:numel(its)
    it = (its(k));
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
% x = 1:10000;
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
plot(x,y);
hold on
plot(x,yf,'.r');
title(['Iteration number',num2str(its(k))]);
legend('Data',' Neural Network');
hold off
drawnow
M(k) = getframe(gcf);
end