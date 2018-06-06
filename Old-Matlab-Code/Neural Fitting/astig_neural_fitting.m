%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Loraxian Transform Preliminary
%
% A transform to take a pixelated image of a psf into the localized
% coordinates
%
% AJN 12/21/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; 
clearvars; clc;

%% Recreate functional fitting
% load('3d.mat','x','sigx','sigy')

xsub = (0:0.01:4).';
ysub = sin(xsub);
x = xsub;
for i = 1:numel(ysub(1,:))
plot(xsub,ysub(:,i))
hold on
end

%% Neural Network Stuff 
% Layer information
inputs = numel(xsub(1,:));
outputs = numel(ysub(1,:));
hiddens = 30;
it = 100;

% Fitting Stuff
lambda = 5;
trainper = 0.6;
epsint = 6;

% build random thetas
% theta1 = rand(hiddens, inputs + 1)*2*epsint - epsint;
% theta2 = rand(outputs, hiddens +1)*2*epsint - epsint;
theta1 = ones(hiddens, inputs +1);
theta2 = ones(outputs, hiddens +1);

% Feature scaling
% for i = 1:numel(ysub(1,:))
% ymi(i) = min(ysub(:,i));
% ysc(i) = (max(ysub(:,i)) - min(ysub(:,i)));
% ys(:,i) = (ysub(:,i) - ymi(i))/ysc(i);
% end
ys = ysub;
disp('here')
[theta1, theta2] = func_lorax_neural_teach(theta1, theta2, xsub, ys, lambda, it, trainper);

% trained neural calc
% x = 1:10000;
x = x(:);
X1 = [ ones(numel(x(:,1)),1), x];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
a2 = [ones(numel(z2(:,1)),1), relu(z2)];
z3 = a2*theta2.';
% a3 = sigmoid(z3);
a3 = relu(z3);

% unscale
for i = 1:numel(ysub(1,:))
% yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
yf(:,i) = a3(:,i) ;
end

hold on
plot(x,yf,'.r');
legend('Data',' Neural Network');
hold off