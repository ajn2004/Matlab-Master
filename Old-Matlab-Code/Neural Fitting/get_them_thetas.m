function [theta1, theta2,ysc,ymi] = get_them_thetas(x,y,hids,l,epsint)
% x and y must be inserted as coloumn vectors
xsub = x;
ysub = y;

%% Neural Network Stuff 
% Layer information
inputs = numel(xsub(1,:));
outputs = numel(ysub(1,:));
hiddens = hids;
it = 10000;

% Fitting Stuff
lambda = l;
trainper = 1;
% epsint = 0.12;

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

% unscale
% for i = 1:numel(ysub(1,:))
% % yf(:,i) = a3(:,i) * ysc(i) + ymi(i);
% end
