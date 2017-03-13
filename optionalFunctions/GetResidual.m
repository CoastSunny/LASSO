function res = GetResidual(Y,X,betaold)
res = Y;
G = numel(betaold);
for i = 1:G
    res = res - X{i}*betaold{i};
end