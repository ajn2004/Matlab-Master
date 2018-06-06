function a3 = func_eval_NN(x,t1,t2)
X1 = [ ones(numel(x(:,1)),1), x];
z2 = X1*t1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*t2.';
a3 = sigmoid(z3);