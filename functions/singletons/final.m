clear; clc;
%addpath /home/mhl2111/MATLAB:/home/mhl2111/MATLAB/PCA:/home/mhl2111/MATLAB/Benoit:/home/mhl2111/MATLAB/forwardData/ROI_correlation_many_subjects:/home/mhl2111/MATLAB/forwardData:/home/mhl2111/MATLAB/Diagnostics
addpath(genpath('/Volumes/Denali_4D2/kohler/LASSO/'));    
load Subject_48_initialization
numComponents = 5;
G = 18;
RUNS = 1;
TOL = 1e-5;
MAX_ITER = 1e5;
VAR_EXPLAINED = 0.99;

[Xlist, Vlist, grpSize] = get_X_list(numComponents);

subjectsG = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];

auc = zeros(numel(subjectsG), RUNS);
relativeEnergy = zeros(numel(subjectsG), RUNS);
for numSubjects = 25:numel(subjectsG)
    fprintf('Working on %d subjects\n', numSubjects);
    rand('seed', 03182012);
    for run = 1:RUNS
        fprintf('Run number: %d\n', run);
        subjects = subjectsG(randsample(numel(subjectsG), numSubjects));
        idx = randsample(G, 2);
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
        %generate Y and signal
        rois = cell(1, numSubjects);
        for N = 1:numSubjects
            load(['forwardAndRois-skeri' num2str(subjects(N))]);
            load(['ROI_correlation_subj_' num2str(subjects(N))]);
            rois{N} = ROIs.ndx;
            [Y(128*(N-1)+1:128*N,:), ~, signal{N}] = GenerateData(ROIs,idx,VertConn,fwdMatrix,0);
        end
        
        %get number of columns of v
        [~, d, v] = svd(Y);
        d = diag(d);
        numCols = 2;
%         while d(1:numCols)'*d(1:numCols) / (d'*d) <= 0.99
%             numCols = numCols + 1;
%         end
        
        %transformed problem
        Ytrans = Y * v(:,1:numCols);
        for i = 1:numSubjects
            signal{i} = signal{i} * v(:, 1:numCols);
        end
        
        %sequence of lambda values
        lambdaMax = 0;
        for i = 1:size(X,2)
            lambdaMax = max(lambdaMax, norm(X(:,i)'*Ytrans));
        end
        lambda = lambdaMax * (0.05.^(0:1/49:1));
        
        %fitting
        betaInit = zeros(size(X,2), numCols);
        beta = cell(1, numel(lambda));
        objValues = cell(1, numel(lambda));
        normY = norm(Ytrans, 'fro')^2;
        for i = 1:numel(lambda)
            [beta{i} objValues{i}] = get_solution(X, Ytrans, betaInit, lambda(i), TOL, MAX_ITER);
            if norm(X*beta{i}, 'fro')^2/normY > VAR_EXPLAINED
                break;
            end
            betaInit = beta{i};
        end
        
        %convert beta back into original space
        betaOriginal = V * beta{i};
        
        aucTemp = zeros(1,numSubjects);
        relativeEnergyTemp = zeros(1,numSubjects);
        grpsizes = zeros(1,G);
        for j = 1:G
            for i = 1:numSubjects
                grpsizes(j) = grpsizes(j) + numel(rois{i}{j});
            end
        end
        for i = 1:numSubjects
            [aucTemp(i) relativeEnergyTemp(i)] = get_metrics(betaOriginal(return_index(grpsizes, rois, i),:), signal{i}, VertConn, rois{i});
        end
        
        auc(numSubjects, run) = mean(aucTemp);
        relativeEnergy(numSubjects, run) = mean(relativeEnergyTemp);
        
    end
end

%save('results.mat', 'auc', 'relativeEnergy');
%exit;
