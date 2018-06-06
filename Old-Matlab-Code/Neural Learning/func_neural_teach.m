function [ Theta1 Theta2 ] = func_neural_teach( int_Theta1, int_Theta2, X, Y, lambda, it, trainper )
% Neural Teach is a function that learns new theta values for x and y data
int_thetas = [int_Theta1(:); int_Theta2(:)];
hidden_layer = numel(int_Theta1(:,1));
input_layer = numel(int_Theta1(1,:))-1;
output = numel(int_Theta2(:,1));

split = round(trainper * numel(X(:,1)));
% conv = 0.001; %convergence factor, ask J to be within 100*conv % of itself between iterations before breaking analysis
% [fin_thetas, i, Jout] = newtraphs(int_thetas, X, y.', lambda, conv); %apply netwon raphson method to determine proper values for theta parameters
randind = randperm(numel(X(:,1)));
Xtrain = X(randind(1:split),:);
Ytrain = Y(randind(1:split));
options = optimset('MaxIter', 100);
[thetas] = fmincg(@(t)(neuralcostfunc(t, input_layer, hidden_layer, output, Xtrain, Ytrain, lambda)), int_thetas, options);

Theta1 = reshape(thetas(1: hidden_layer*(input_layer +1)), hidden_layer, input_layer + 1);
Theta2 = reshape(thetas(hidden_layer*(input_layer +1) +1: end), output, hidden_layer + 1);

%% Determine accuracy
X1 = [ ones(numel(X(:,1)),1), X];
z2 = X1*Theta1.';
a2 = [ones(numel(z2(:,1)),1), sigmoid(z2)];
z3 = a2*Theta2.';
a3 = sigmoid(z3);
Xtest = X(randind(split+1:end),:);
Ytest = Y(randind(split+1:end));
ptrain = mean(a3 > 0.5 == Y) %
xt1 = [ ones(numel(Xtest(:,1)),1), Xtest];
zt2 = xt1*Theta1.';
at2 =[ones(numel(zt2(:,1)),1), sigmoid(zt2)];
zt3 = at2*Theta2.';
at3 = sigmoid(zt3);
ptest = mean(at3>0.5 == Ytest)

% save([current_path, 'neural_theta.mat'],'Theta1','Theta2');
if it < 10
    save(['Neural_thetas_it_0',num2str(it),'.mat'],'Theta1','Theta2');
else
    save(['Neural_thetas_it_',num2str(it),'.mat'],'Theta1','Theta2');
end
end

