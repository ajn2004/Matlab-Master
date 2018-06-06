function [J, grad] = neural_multi_class_cost(thetas,i,h,o, x, y, index, lambda)

% Reshape thetas
Theta1 = reshape( thetas(1:(i+1)*h) , h,i+1);
Theta2 = reshape(thetas((i+1)*h+1:end), o, h + 1);
m = numel(x(:,1)); % number of examples
n = numel(Theta2(:,1));
Y = zeros(m,n);
Y(:,end) = 1; % last component is a "trash" component

for i = 1:n-1 % convert vectors into 0 and 1s
    
    ind = y == index(i);
    Y(ind,i) = 1;
    Y(ind,end) = 0;

end

% Preallocate J
J = 0;
T1_grad = zeros(size(Theta1));
T2_grad = zeros(size(Theta2));

% Feed Forward and get current values for activation of the second layer
X = [ones(m,1), x];
z1 = X*Theta1.';
a1 = sigmoid(z1);
a1f = [ones(m,1), a1];
z2 = a1f*Theta2.';
a2 = sigmoid(z2);

% calculate the cost of all of these
J = (1/m)*sum(sum( -Y.* log(a2) - (1-Y).*log(1 - a2))) + lambda/(2*m)*(sum(sum(Theta1(:,2:end).^2)) + sum(sum(Theta2(:,2:end).^2)));

% Backpropogation
del3 = (a2 - Y);
d2 = del3*Theta2(:,2:end);
del2 = d2.*sigmoidGradient(z1);

% Given del 2 we can calculate the gradient for neural parameters.
T1_grad = (del2.'*X)/m;
T1_grad(:,2:end) = T1_grad(:,2:end) + lambda/m*Theta1(:,2:end);
T2_grad = (del3.'* a1f)/m;
T2_grad(:,2:end) = T2_grad(:,2:end) + lambda/m*Theta2(:,2:end);
grad = [T1_grad(:) ; T2_grad(:)];

end