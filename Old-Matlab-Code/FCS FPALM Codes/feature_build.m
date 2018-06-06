function x = feature_build(X, i)

x = ones(numel(X(:,1)),i+1);

for j = 2:i+1
    x(:,j) = X.^(j-1);

end