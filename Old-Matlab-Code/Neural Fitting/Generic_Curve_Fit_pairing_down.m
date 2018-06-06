%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic Curve Fitting
%
% A transform to take a pixelated image of a psf into the localized
% coordinates
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc;
% 
% %% Use this section to create your data set 
% % load('composit.mat','x','drifts')
% 
% %% 
% xsub = (-5:0.01:5).';
% ysub = exp(-xsub.^2/(2*2^2));
close all
%% Plot initial data
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
lambda = 0.001;
trainper = 1;
epsint = 0.82;

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
x = xsub;
a3 = func_eval_NN(x,theta1,theta2);


% unscale
for i = 1:numel(ysub(1,:))
yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
end

hold on
plot(x,yf,'.r');
legend('Data',' Neural Network');