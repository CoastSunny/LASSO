function [result, v] = get_principal_components(X, numComponents)

if isempty(X)
    result = NaN;
    v = NaN;
else
    X = scal(X, mean(X));

    [~, ~, v] = svd(X);

    v = v(:, 1:numComponents);
    result = X*v;
end

