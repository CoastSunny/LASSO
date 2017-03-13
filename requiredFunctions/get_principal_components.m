function [result, v] = get_principal_components(X, numComponents)

%scal not a matlab builtin and I don't know where its from
%X = scal(X, mean(X));
%fixed to use the bsxfun which is a matlab builtin for a few years now.
X = bsxfun(@minus,X,mean(X));

[~, ~, v] = svd(X);

v = v(:, 1:numComponents);
result = X*v;

