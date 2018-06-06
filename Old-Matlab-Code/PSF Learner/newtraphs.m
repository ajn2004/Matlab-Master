function [theta_out, i, outJ] = newtraphs(theta, x, y, lambda, conv)
% this function will iteritavely apply the Newton-Raphson method to
% minimize the cost function

theta_out = theta.*0;
m = numel(y);
i = 2;
conv = abs(conv);
J = [2, 20];
while i == 1 || J(i) > 0.01 || (J(i) > (1+conv)*J(i-1) || J(i) < (1-conv)*J(i-1))
    i = i+1;
    djdt = m^-1 * x.'*(sigmoid(x*theta) - y);
    temp = lambda* theta / m;
    temp(1) = 0;
    djdt = djdt + temp;
    
    d2jdt2 = m^-1 *x.^2.' * sigmoid(x*theta).^2*sum(exp(-x*theta));
    temp = temp*0 + lambda / m;
    temp(1) = 0;
    d2jdt2 = d2jdt2 + temp;
    
    theta = theta - djdt./d2jdt2;
    [J(i) grad] = costfunc(theta, x, y, lambda);
    
end
theta_out = theta(:);
outJ = J(i);
end