function [J, Grads] = neuralcostfunc(thetas, ...
                                     in_layer, ...
                                     hidden_layer, ...
                                     output, ...
                                     X, y, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Cost Function
% 
% This is the function which will determine the cost of parameters  in thetas
% Thetas will be a 1 dimensional vector that has to be reshaped inside the
% function, in layer tells us how  X is the raw data and
% y is the indicator of the data. Lambda is the regulation parameters to
% help prevent overfitting
%
% AJN 10/7/15
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% reshape theta1 and theta2
Theta1 = reshape(thetas(1:hidden_layer*(in_layer+1)), hidden_layer, in_layer+1);
Theta2 = reshape(thetas(hidden_layer*(in_layer+1)+1:end), output, hidden_layer+1);

m = numel(X(:,1));

% initialize outputs
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));
J = 0;

%Build values for cost and gradients
x = [ones(numel(X(:,1)),1), X];
z2 = x*Theta1.';
a2 = [ones(numel(z2(:,1)),1),sigmoid(z2)];
z3 = a2*Theta2.';
a3 = sigmoid(z3);

% Regularized cost function
J = (1/m)*sum(sum( -y.* log(a3) - (1-y).*log(1 - a3))) + lambda/(2*m)*(sum(sum(Theta1(:,2:end).^2)) + sum(sum(Theta2(:,2:end).^2)));

% Gradient of thetas
del3 = a3 - y ;
d2 = del3*Theta2(:,2:end);
del2 = d2.*sigmoidGradient(z2);

 Theta1_grad = (del2.'*x)/m;
 Theta1_grad(:,2:end) = Theta1_grad(:,2:end) + lambda/m*Theta1(:,2:end);
 Theta2_grad = (del3.'* a2)/m;
 Theta2_grad(:,2:end) = Theta2_grad(:,2:end) + lambda/m*Theta2(:,2:end);

 % roll up thetas for export
Grads = [Theta1_grad(:); Theta2_grad(:)];