function y = sigmoid(x)
% this function calculates the value of the sigmoid for x

y = (1+ exp(-x)).^-1;
end