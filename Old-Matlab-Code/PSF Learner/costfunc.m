function [J, grad] = costfunc(theta, x, y, lambda)
% this function will compute the cost of fitting x with parameters theta

J = 0;
grad = zeros(numel(theta),1);
m = numel(y(:,1)); % number of training examples


J = m^-1 * sum(- y.* log(sigmoid(x*theta)) - (1-y).*log(1-sigmoid(x*theta)+.000000000000000000001)) + lambda / (2*m) * sum(theta(2:end).^2);

grad = m^-1 * x.'*(sigmoid(x*theta) - y);
temp = lambda/m * theta;
temp(1) = 0;
grad = grad + temp;
end