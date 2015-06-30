clear; clc;
%addpath /home/mhl2111/MATLAB:/home/mhl2111/MATLAB/PCA/frobenius:/home/mhl2111/MATLAB/Benoit:/home/mhl2111/MATLAB/forwardData/ROI_correlation_many_subjects:/home/mhl2111/MATLAB/forwardData:/home/mhl2111/MATLAB/Diagnostics
load Subject_48_initialization
numComponents = 5;
numCols = 2;
nLambda = 50;
nLambdaRidge = 50;
alpha = 1.0817e4;
G = 18;
RUNS = 1;
MAX_ITER = 1e6;
SNR = 0.1;
numSubjects = 5;
idx = [3, 12];
rng(29012013, 'twister');

%Xlist is the list of derived forward matrices, each of dimension (n x numComponents). Vlist contains the matrices used to reverse the PCA transformation
[Xlist, Vlist] = get_X_list(numComponents);

%subject IDs, which correspond to the file names
subjectsG = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];
totalSubjects = numel(subjectsG);

aucCloseRidge = zeros(1, nLambda);
aucFarRidge = zeros(1, nLambda);
mseRidge = zeros(1, nLambda);
energyRidge = zeros(1, nLambda);
rssRidgeRuns = zeros(1, nLambda);
aucClose = zeros(1, nLambda);
aucFar = zeros(1, nLambda);
mse = zeros(1, nLambda);
energy = zeros(1, nLambda);
rss = zeros(1, nLambda);

subjects = subjectsG(randsample(numel(subjectsG), numSubjects)); %sample subjects at random
signal = cell(1, numSubjects);
Y = zeros(numSubjects*128, 91);
%generate X and V
X = [];
V = [];
for g = 1:G
    tempX = [];
    tempV = [];
    for N = 1:numSubjects
        tempX = blkdiag(tempX, Xlist{subjectsG==subjects(N)}{g});
        tempV = blkdiag(tempV, Vlist{subjectsG==subjects(N)}{g}(:,1:numComponents));
    end
    X = [X, tempX];
    V = blkdiag(V, tempV);
end
grpSizes = numComponents*numSubjects*ones(1,G);
indices = get_indices(grpSizes);
penalties = get_group_penalties(X, indices);
%generate Y and signal
rois = cell(1, numSubjects);    
stackedForwards = [];
phase = randi([2,10], 1);
for N = 1:numSubjects
    load(['forwardAndRois-skeri' num2str(subjects(N))]);
    load(['ROI_correlation_subj_' num2str(subjects(N))]);
    rois{N} = ROIs.ndx;
    [Y(128*(N-1)+1:128*N,:), ~, signal{N}, noise] = GenerateData(ROIs,idx,VertConn,fwdMatrix,SNR, phase);
    stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[ROIs.ndx{:}]));
end
%done with generating Y and X

%center Y, X, stackedForwards
X = scal(X, mean(X));
stackedForwards = scal(stackedForwards, mean(stackedForwards));
Y = scal(Y, mean(Y));
n = numel(Y);
ssTotal = norm(Y, 'fro')^2 / n;

%use first 2 columns of v as time basis
[~, ~, v] = svd(Y);

%transformed problem
Ytrans = Y * v(:, 1:numCols);
Ytrans = scal(Ytrans, mean(Ytrans));

%minumum norm solution
[betaRidgeBest, betaRidge, lambdaMinNorm, gcvErrorRidge, lambdaGridRidge] = minimum_norm(stackedForwards, Y, nLambdaRidge);
lambdaMinNorm = lambdaMinNorm^2;
lambdaGridRidge = lambdaGridRidge.^2;
rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaRidgeBest, 'fro')^2/n) / ssTotal;
rsquaredRidge = zeros(1, nLambdaRidge);
for i = 1:nLambdaRidge
    rsquaredRidge(i) = 1 - (norm(Y-stackedForwards*betaRidge{i}, 'fro')^2/n) / ssTotal;
end

%sequence of lambda values
lambdaMax = 0;
for i = 1:G
    lambdaMax = max(lambdaMax, norm(X(:,indices{i})'*Ytrans, 'fro')/penalties(i));
end
lambdaMax = lambdaMax + 1e-4;
lambdaGrid = lambdaMax * (0.001.^(0:1/(nLambda-1):1));
tol = min(penalties) * lambdaGrid(end) * 1e-5;
if alpha > 0
    tol = min([tol, 2*alpha*1e-5]);
end

ridgeSizes = zeros(1, numSubjects);
for i = 1:numSubjects
    ridgeSizes(i) = numel([rois{i}{:}]);
end                    
ridgeRange = [0 cumsum(ridgeSizes)];

roiSizes = zeros(1,G); %total size of each region across all subjects
for j = 1:G
    for i = 1:numSubjects
        roiSizes(j) = roiSizes(j) + numel(rois{i}{j});
    end
end

%ols fit
betaOls = (X'*X + alpha*eye(size(X,2))) \ (X'*Ytrans);

%fitting
fprintf('Working on %d subjects\n', numSubjects);
betaInit = zeros(size(X,2), numCols);
beta = cell(1, nLambda);
objValues = cell(1, nLambda);
gcvError = zeros(1, nLambda);
df = zeros(1, nLambda);
rsquared = zeros(1, nLambda);
converged = zeros(1, nLambda);
indexer = [];
for i = 1:numSubjects
    indexer = [indexer, return_index(roiSizes, rois, i)];
end
fprintf('i\tlambda\tobj      \tniter\tgcvError\tdf\tconverged\n');
fprintf('--------------------------------------------------------------------------------\n');
for i = 1:nLambda
    [beta{i}, objValues{i}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(i), alpha, tol, MAX_ITER, penalties, indices);
    converged(i) = check_solution(X, beta{i}, res, indices, penalties, lambdaGrid(i), alpha, tol);
    betaInit = beta{i};
    beta{i} = V * beta{i} * v(:,1:numCols)'; %transform back to original space (permuted forward matrices)
    rss(i) = norm(Y-stackedForwards*beta{i}(indexer, :), 'fro')^2 / n;
    rsquared(i) = 1 - rss(i) / ssTotal;
    [gcvError(i), df(i)] = compute_gcv(rss(i), betaInit, betaOls, grpSizes, n);
    fprintf('%d\t%2g\t%2g\t%d\t%2g\t%2g\t%5d\n', i, lambdaGrid(i), objValues{i}(end), numel(objValues{i}), gcvError(i), df(i), converged(i));
end
[~, bestIndex] = min(gcvError);

%compute auc, mse for ridge fit
aucCloseTemp = zeros(1, numSubjects);
aucFarTemp = zeros(1, numSubjects);
mseTemp = zeros(1, numSubjects);
energyTemp = zeros(1, numSubjects);
for j = 1:nLambdaRidge
    for i = 1:numSubjects
        [~, aucCloseTemp(i), aucFarTemp(i), mseTemp(i), ~, energyTemp(i)] = get_metrics(betaRidge{j}(ridgeRange(i)+1:ridgeRange(i+1), :), signal{i}, VertConn, rois{i});
    end
    aucCloseRidge(j) = mean(aucCloseTemp);
    aucFarRidge(j) = mean(aucFarTemp);
    mseRidge(j) = mean(mseTemp);
    energyRidge(j) = mean(energyTemp);
end

%compute auc, mse for lasso
for j = 1:nLambda
    for i = 1:numSubjects
        [~, aucCloseTemp(i), aucFarTemp(i), mseTemp(i), ~, energyTemp(i)] = get_metrics(beta{j}(return_index(roiSizes, rois, i), :), signal{i}, VertConn, rois{i});
    end
    aucClose(j) = mean(aucCloseTemp);
    aucFar(j) = mean(aucFarTemp);
    mse(j) = mean(mseTemp);
    energy(j) = mean(energyTemp);
end

%compute average metrics
T = size(Y, 2);
truth = zeros(1, G);
truth(idx) = 1;
regionActivityRidge = zeros(G, T);
regionActivity = zeros(G, T);
for s = 1:numSubjects
    range = zeros(1, G);
    for i = 1:G
        range(i) = numel(rois{s}{i}); %correct indexing of rois for each subject s
    end
    range = [0 cumsum(range)];
    tempRidge = betaRidgeBest(ridgeRange(s)+1:ridgeRange(s+1), :);
    temp = beta{bestIndex}(return_index(roiSizes, rois, s), :);
    for i = 1:G
        regionActivityRidge(i, :) = regionActivityRidge(i, :) + mean(tempRidge(range(i)+1:range(i+1), :));
        regionActivity(i, :) = regionActivity(i, :) + mean(temp(range(i)+1:range(i+1), :));
    end
end
regionActivityRidge = regionActivityRidge / numSubjects;
regionActivity = regionActivity / numSubjects;
temp = zeros(1, T);
for i = 1:T
    temp(i) = temp(i) + rocArea(abs(regionActivityRidge(:, i)), truth);
end
aucRidgeAvg = mean(temp);
temp = zeros(1, T);
for i = 1:T
    temp(i) = temp(i) + rocArea(abs(regionActivity(:, i)), truth);
end
aucAvg = mean(temp);

%compute least squares fit on reduced problem (end of path solution)
betaLs = V*((X'*X+alpha*eye(size(X,2)))\(X'*Ytrans))*v(:,1:numCols)';
betaLs = betaLs(indexer, :);
rsquaredOls = 1 - (norm(Y-stackedForwards*betaLs, 'fro')^2/n)/ssTotal;
dfOls = 2 + G*numComponents*numSubjects*numCols; %use same df formula for ols fit, not the trace
gcvErrorOls = (norm(Y-stackedForwards*betaLs, 'fro')^2/n) / (1-dfOls/n)^2;
for i = 1:numSubjects
    [~, aucCloseTemp(i), aucFarTemp(i), mseTemp(i), ~, energyTemp(i)] = get_metrics(betaLs(ridgeRange(i)+1:ridgeRange(i+1), :), signal{i}, VertConn, rois{i});
end
aucCloseOls = mean(aucCloseTemp);
aucFarOls = mean(aucFarTemp);
mseOls = mean(mseTemp);
energyOls = mean(energyTemp);

aucClose = [aucClose, aucCloseOls];
aucFar = [aucFar, aucFarOls];
mse = [mse, mseOls];
energy = [energy, energyOls];
rsquared = [rsquared, rsquaredOls];
gcvError = [gcvError, gcvErrorOls];

%check variance of least squares fit
variance = zeros(1, nLambdaRidge);
for i = 1:nLambdaRidge
    temp = V*((X'*X+lambdaGridRidge(i)*eye(size(X, 2)))\X');
    variance(i) = norm(temp, 'fro')^2;
end

file = ['~/Dropbox/michael_research/eeg/writeup/graphics/', num2str(numSubjects), '_subject.mat'];
save(file, 'aucClose', 'aucCloseRidge', 'aucFar', 'aucFarRidge', 'mse', 'mseRidge', 'energy', 'energyRidge', 'aucAvg', 'aucRidgeAvg', 'rsquared', 'rsquaredMinNorm', 'rsquaredRidge', 'lambdaGridRidge', 'variance', 'converged', 'gcvError', 'gcvErrorRidge', 'df');