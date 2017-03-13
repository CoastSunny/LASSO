clear; clc;
%addpath /home/mhl2111/MATLAB:/home/mhl2111/MATLAB/PCA/frobenius:/home/mhl2111/MATLAB/Benoit:/home/mhl2111/MATLAB/forwardData/ROI_correlation_many_subjects:/home/mhl2111/MATLAB/forwardData:/home/mhl2111/MATLAB/Diagnostics
load Subject_48_initialization
G = 18;
RUNS = 100;
SNR = 0.1;

mseGap = zeros(8, RUNS);
gap = zeros(8, RUNS);

for numComponents = 1:8

    %Xlist is the list of derived forward matrices, each of dimension (n x numComponents). Vlist contains the matrices used to reverse the PCA transformation
    [Xlist, Vlist] = get_X_list(numComponents);

    %subject IDs, which correspond to the file names
    subjectsG = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];
    totalSubjects = numel(subjectsG);

    for numSubjects = 1%1:totalSubjects
        rng(29012013, 'twister');
        for run = 1:RUNS
            %idx = randsample(G, 2); %sample regions at random
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
            %normalize X to be orthonormal within each group
            grpSizes = numComponents*numSubjects*ones(1,G);
            indices = get_indices(grpSizes);
            normalizer = zeros(1, numComponents*numSubjects);
            for i = 1:G
                normalizer(indices{i}) = 1./sqrt(diag(X(:,indices{i})'*X(:,indices{i})));
            end
            X = X * diag(normalizer);

            %generate Y and signal
            rois = cell(1, numSubjects);    
            stackedForwards = [];
            phase = randi([2,10], 1);
            for N = 1:numSubjects
                load(['forwardAndRois-skeri' num2str(subjects(N))]);
                load(['ROI_correlation_subj_' num2str(subjects(N))]);
                rois{N} = ROIs.ndx;
                [Y(128*(N-1)+1:128*N,:), ~, signal{N}, noise] = GenerateData(ROIs,idx,VertConn,fwdMatrix,SNR, phase);
            end
            
            a = signal{1}([rois{1}{:}], :);
            mseGap(numComponents, run) = norm(a - V*V'*a, 'fro')^2 / numel(signal{1});
            
            [~, ~, v] = svd(Y);
            gap(numComponents, run) = norm(signal{1} - signal{1}*v(:,1:numComponents)*v(:,1:numComponents)', 'fro')^2 / numel(signal{1});
            
        end
    end
end
    
mseGapMean = mean(mseGap, 2);
mseGapSd = std(mseGap, 1, 2);
plot(1:8, mseGapMean, 1:8, mseGapMean-mseGapSd, '--b', 1:8, mseGapMean+mseGapSd, '--b', 'LineWidth', 2);
ylabel('$\|\beta-PP^T\beta\|_F^2$', 'interpreter', 'latex');
xlabel('Number of principal components');

figure;
gapMean = mean(gap, 2);
gapSd = std(gap, 1, 2);
plot(1:8, gapMean, 1:8, gapMean-gapSd, '--b', 1:8, gapMean+gapSd, '--b', 'LineWidth', 2);
ylabel('$\|\beta-\beta VV^T\|_F^2$', 'interpreter', 'latex');
xlabel('Number of right singular vectors of Y');
