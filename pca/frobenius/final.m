clear; clc;
%addpath /home/mhl2111/MATLAB:/home/mhl2111/MATLAB/PCA/frobenius:/home/mhl2111/MATLAB/Benoit:/home/mhl2111/MATLAB/forwardData/ROI_correlation_many_subjects:/home/mhl2111/MATLAB/forwardData:/home/mhl2111/MATLAB/Diagnostics
addpath(genpath('/Volumes/Denali_4D2/kohler/LASSO'));
load Subject_48_initialization
numComponents = 5;
numCols = 2;
alpha = 1.0817e4;
nLambda = 30;
nLambdaRidge = 50;
G = 18;
RUNS = 20;
MAX_ITER = 1e6;
SNR = 0.1;

%Xlist is the list of derived forward matrices, each of dimension (n x numComponents). Vlist contains the matrices used to reverse the PCA transformation
[Xlist, Vlist] = get_X_list(numComponents);

%subject IDs, which correspond to the file names
subjectsG = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];
totalSubjects = numel(subjectsG);

%Loop through number of subjects. Generate random signal for regions using a randomly generated phase that is the same across all subjects
aucCloseMinNormRuns = zeros(RUNS, totalSubjects);
aucFarMinNormRuns = zeros(RUNS, totalSubjects);
mseMinNormRuns = zeros(RUNS, totalSubjects);
energyMinNormRuns = zeros(RUNS, totalSubjects);
aucCloseRuns = zeros(RUNS, totalSubjects);
aucFarRuns = zeros(RUNS, totalSubjects);
mseRuns = zeros(RUNS, totalSubjects);
energyRuns = zeros(RUNS, totalSubjects);

aucMinNormRunsAvg = zeros(RUNS, totalSubjects);
aucRunsAvg = zeros(RUNS, totalSubjects);
%matlabpool 5;
for numSubjects = 25;%[1, 2, 4, 8, 16, 25]
    fprintf('Working on %d subjects\n', numSubjects);
    rng(29012013, 'twister');
    for run = 1:RUNS
        fprintf('Run number: %d\n', run);
        idx = [3, 12];
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
            structure = load(['forwardAndRois-skeri' num2str(subjects(N))]);
            fwdMatrix = structure.fwdMatrix;
            structure = load(['ROI_correlation_subj_' num2str(subjects(N))]);
            ROIs = structure.ROIs;
            rois{N} = ROIs.ndx;
            [Y(128*(N-1)+1:128*N,:), ~, signal{N}, noise] = GenerateData(ROIs,idx,VertConn,fwdMatrix,SNR, phase);
            stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[ROIs.ndx{:}]));
        end
        %done with generating Y and X
        
        %center Y, X, stackedForwards
        Y = scal(Y, mean(Y));
        X = scal(X, mean(X));
        stackedForwards = scal(stackedForwards, mean(stackedForwards));
        n = numel(Y);
        ssTotal = norm(Y, 'fro')^2 / n;
        
        %use first 2 columns of v as time basis
        [~, ~, v] = svd(Y);
        Ytrans = Y * v(:, 1:numCols);
        Ytrans = scal(Ytrans, mean(Ytrans));
        
        %minumum norm solution
        [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Y, nLambdaRidge);
        lambdaMinNorm = lambdaMinNorm^2;
        rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotal;
        
        %sequence of lambda values
        lambdaMax = 0;
        for i = 1:G
            lambdaMax = max(lambdaMax, norm(X(:,indices{i})'*Ytrans, 'fro')/penalties(i));
        end
        lambdaMax = lambdaMax + 1e-4;
        lambdaGrid = lambdaMax * (0.01.^(0:1/(nLambda-1):1));
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
        betaInit = zeros(size(X,2), numCols);
        beta = cell(1, nLambda);
        objValues = cell(1, nLambda);
        gcvError = zeros(1, nLambda);
        df = zeros(1, nLambda);
        rsquared = zeros(1, nLambda);
        indexer = [];
        for i = 1:numSubjects
            indexer = [indexer, return_index(roiSizes, rois, i)];
        end
        for i = 1:nLambda
            [beta{i}, objValues{i}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(i), alpha, tol, MAX_ITER, penalties, indices);
            betaInit = beta{i};
            beta{i} = V * beta{i} * v(:,1:numCols)'; %transform back to original space (permuted forward matrices)
            rss = norm(Y-stackedForwards*beta{i}(indexer, :), 'fro')^2 / n;
            [gcvError(i), df(i)] = compute_gcv(rss, betaInit, betaOls, grpSizes, n);
        end
        [~, bestIndex] = min(gcvError);
            
        %compute auc, mse for min norm
        aucCloseTemp = zeros(1, numSubjects);
        aucFarTemp = zeros(1, numSubjects);
        mseTemp = zeros(1, numSubjects);
        energyTemp = zeros(1, numSubjects);
        for i = 1:numSubjects
            [~, aucCloseTemp(i), aucFarTemp(i), mseTemp(i), ~, energyTemp(i)] = get_metrics(betaMinNorm(ridgeRange(i)+1:ridgeRange(i+1), :), signal{i}, VertConn, rois{i});
        end
        aucCloseMinNormRuns(run, numSubjects) = mean(aucCloseTemp);
        aucFarMinNormRuns(run, numSubjects) = mean(aucFarTemp);
        mseMinNormRuns(run, numSubjects) = mean(mseTemp);
        energyMinNormRuns(run, numSubjects) = mean(energyTemp);
        
        %compute auc, mse for lasso using optimal lambda
        for i = 1:numSubjects
            [~, aucCloseTemp(i), aucFarTemp(i), mseTemp(i), ~, energyTemp(i)] = get_metrics(beta{bestIndex}(return_index(roiSizes, rois, i), :), signal{i}, VertConn, rois{i});
        end
        aucCloseRuns(run, numSubjects) = mean(aucCloseTemp);
        aucFarRuns(run, numSubjects) = mean(aucFarTemp);
        mseRuns(run, numSubjects) = mean(mseTemp);
        energyRuns(run, numSubjects) = mean(energyTemp);
        
        %compute average of average metrics
        T = size(Y, 2);
        truth = zeros(1, G);
        truth(idx) = 1;
        regionActivityMinNorm = zeros(G, T);
        regionActivity = zeros(G, T);
        for s = 1:numSubjects
            range = zeros(1, G);
            for i = 1:G
                range(i) = numel(rois{s}{i}); %correct indexing of rois for each subject s
            end
            range = [0 cumsum(range)];
            tempMinNorm = betaMinNorm(ridgeRange(s)+1:ridgeRange(s+1), :);
            temp = beta{bestIndex}(return_index(roiSizes, rois, s), :);
            for i = 1:G
                regionActivityMinNorm(i, :) = regionActivityMinNorm(i, :) + mean(tempMinNorm(range(i)+1:range(i+1), :));
                regionActivity(i, :) = regionActivity(i, :) + mean(temp(range(i)+1:range(i+1), :));
            end
        end
        regionActivityMinNorm = regionActivityMinNorm / numSubjects;
        regionActivity = regionActivity / numSubjects;
        temp = zeros(1, T);
        for i = 1:T
            temp(i) = temp(i) + rocArea(abs(regionActivityMinNorm(:, i)), truth);
        end
        aucMinNormRunsAvg(run, numSubjects) = mean(temp);
        temp = zeros(1, T);
        for i = 1:T
            temp(i) = temp(i) + rocArea(abs(regionActivity(:, i)), truth);
        end
        aucRunsAvg(run, numSubjects) = mean(temp);
    end
end
%matlabpool close;

save('/Volumes/Denali_4D2/kohler/LASSO/writeup/graphics/pk_simulation_results.mat', 'aucCloseRuns', 'aucFarRuns', 'mseRuns', 'energyRuns','aucCloseMinNormRuns', 'aucFarMinNormRuns', 'mseMinNormRuns', 'energyMinNormRuns', 'aucRunsAvg', 'aucMinNormRunsAvg');