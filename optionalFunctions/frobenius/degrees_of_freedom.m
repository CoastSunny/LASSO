clear; clc
rng(0, 'twister');
load ./pca/frobenius/df_data

runs = 1000;
nLambda = 30;
sigma = 5;
alpha = 0;
[n, p] = size(X);
T = 2;

%generate Y, make first column's variance larger
beta = [normrnd(0, sigma, 10, T); zeros(p-10, T)];
Y = autosc(X*beta);
Y(:, 1) = 5*Y(:, 1);

%compute group penalties
penalties = get_group_penalties(X, indices);
lambdaMax = 0;
for i = 1:18
    lambdaMax = max(lambdaMax, norm(X(:,indices{i})'*Y, 'fro')/penalties(i));
end
lambdaMax = 1.5*lambdaMax;
lambdaGrid = lambdaMax * (0.005.^(0:1/(nLambda-1):1));
tol = min(penalties) * lambdaGrid(end) * 1e-5;
if alpha > 0
    tol = min([tol, 2*alpha*1e-5]);
end

betaNorms = zeros(runs, nLambda, 18);
Ynoisy = cell(1, runs);
df = zeros(runs, nLambda);
dfLs = zeros(1, runs);
matlabpool 5;
parfor run = 1:runs
    if mod(run, 20) == 0
        fprintf('Run number %d\n', run);
    end
    
    %generate noise
    noise = normrnd(0, sigma, n, T);
    Ynoisy{run} = Y + noise;
    Ynoisy{run} = scal(Ynoisy{run}, mean(Ynoisy{run}));
    
    %least squares solution
    YhatLs = X*((X'*X)\(X'*Ynoisy{run}));
    dfLs(run) = trace(scal(YhatLs, mean(YhatLs))'*noise) / sigma^2;
      
    %lasso fit
    betaInit = zeros(p, T);
    for i = 1:nLambda
        [beta, objValues, res] = get_solution_frobenius(X, Ynoisy{run}, betaInit, lambdaGrid(i), alpha, tol, 1e6, penalties, indices);
        betaInit = beta;
        %compute lasso df
        for j = 1:18
            betaNorms(run, i, j) = norm(beta(indices{j},:), 'fro');
        end
        %compute true df
        Yhat = X*beta;
        df(run, i) = trace(scal(Yhat, mean(Yhat))'*noise) / sigma^2;
    end
end
matlabpool close;

save('./pca/frobenius/df_data.mat', 'Y', 'X', 'penalties', 'indices', 'df', 'dfLs', 'Ynoisy', 'betaNorms', 'lambdaGrid', 'tol', 'alpha');

clear;clc
load ./pca/frobenius/df_data
p = 90;
nLambda = 30;
grpsizes = 5*ones(1, 18);
runs = 1000;

%ols solution (optimal value is 1.0817e4)
betaOlsNorms = zeros(runs, 18);
for run = 1:runs
    betaOls = (X'*X + 1.0817e4*eye(p)) \ (X'*Ynoisy{run});
    for i = 1:18
        betaOlsNorms(run, i) = norm(betaOls(indices{i},:), 'fro');
    end
end

%estimated df
dfEstimated = 2*ones(runs, nLambda);
for run = 1:runs
    for i = 1:nLambda
        for j = 1:18
            dfEstimated(run, i) = dfEstimated(run, i) + (betaNorms(run, i, j) > 0) + (2*grpsizes(j)-1)*betaNorms(run, i, j)/betaOlsNorms(run, j);
        end
    end
end

plot(mean(df)+2, mean(dfEstimated), 'x');
hold;
plot([0, max([mean(df)+2, mean(dfEstimated)])], [0, max([mean(df)+2, mean(dfEstimated)])], '--');