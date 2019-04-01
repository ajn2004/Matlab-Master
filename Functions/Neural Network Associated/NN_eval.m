function a3 = NN_eval(x,theta1,theta2)

X1 = [ ones(numel(x(:,1)),1), x];
z2 = X1*theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*theta2.';
a3 = sigmoid(z3);

