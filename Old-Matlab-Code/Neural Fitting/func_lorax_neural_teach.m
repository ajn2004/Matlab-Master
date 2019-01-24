function [ Theta1 Theta2 ] = func_lorax_neural_teach( int_Theta1, int_Theta2, X, Y, lambda, it, trainper )
% Neural Teach is a function that learns new theta values for x and y data
int_thetas = [int_Theta1(:); int_Theta2(:)]; % linearize the thetas
hidden_layer = numel(int_Theta1(:,1));  % Learn how big the hidden layer is
input_layer = numel(int_Theta1(1,:))-1; % Learn how big the input layer is
output = numel(int_Theta2(:,1));        % Learn how big the output layer is

split = round(trainper * numel(X(:,1))); %identifiy how many training examples will be used
% conv = 0.001; %convergence factor, ask J to be within 100*conv % of itself between iterations before breaking analysis
% [fin_thetas, i, Jout] = newtraphs(int_thetas, X, y.', lambda, conv); %apply netwon raphson method to determine proper values for theta parameters
randind = randperm(numel(X(:,1)));
Xtrain = X(randind(1:split),:);
Ytrain = Y(randind(1:split),:);
options = optimset('MaxIter', it);
disp('Training')

% Show Progress

% Adjust theta values to minimize t which is the cost

[thetas] = fmincg(@(t)(neuralcostfunc(t, input_layer, hidden_layer, output, Xtrain, Ytrain, lambda)), int_thetas, options);
disp('Complete')
Theta1 = reshape(thetas(1: hidden_layer*(input_layer +1)), hidden_layer, input_layer + 1);
Theta2 = reshape(thetas(hidden_layer*(input_layer +1) +1: end), output, hidden_layer + 1);
end

