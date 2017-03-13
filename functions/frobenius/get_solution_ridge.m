function [betaRidge, gcvError, df] = get_solution_ridge(Y, X, lambdaGrid)

n = numel(Y);
t = size(Y, 2);
p = size(X, 2);
numLambda = numel(lambdaGrid);
gcvError = zeros(1, numLambda);
df = zeros(1, numLambda);
betaRidge = cell(1, numLambda);
for i = 1:numLambda
    betaRidge{i} = (X'*X+lambdaGrid(i)*eye(p)) \ (X'*Y);
    H = X * ((X'*X+lambdaGrid(i)*eye(p)) \ X');
    rss = norm(Y - H*Y, 'fro')^2 / n;
    df(i) = t * trace(H);
    gcvError(i) = rss / (1-df(i)/n)^2;
end

