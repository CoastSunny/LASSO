function dfEstimated = df_optimizer(shrinkage, structure)

X = structure.X;
Ynoisy = structure.Ynoisy;
indices = structure.indices;
betaNorms = structure.betaNorms;

p = 90;
nLambda = 30;
grpsizes = 5*ones(1, 18);
runs = 1000;

betaOlsNorms = zeros(runs, 18);
for run = 1:runs
    betaOls = (X'*X + shrinkage*eye(p)) \ (X'*Ynoisy{run});
    for i = 1:18
        betaOlsNorms(run, i) = norm(betaOls(indices{i},:), 'fro');
    end
end

dfEstimated = 2*ones(runs, nLambda);
for run = 1:runs
    for i = 1:nLambda
        for j = 1:18
            dfEstimated(run, i) = dfEstimated(run, i) + (betaNorms(run, i, j) > 0) + (2*grpsizes(j)-1)*betaNorms(run, i, j)/betaOlsNorms(run, j);
        end
    end
end

dfEstimated = mean(dfEstimated);