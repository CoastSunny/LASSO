function [result, v] = get_principal_components(X, numComponents)

for i = 1:size(X, 2)
    X(:,i) = X(:,i) - mean(X(:,i));
end

[~, ~, v] = svd(X);

result = X*v(:,1:numComponents);