function [Y, Yhat, YhatMinNorm] = compute_roc(signal, beta, betaMinNorm, rois, grpsizes)

numSubjects = numel(rois);
numLambda = numel(beta);

%get combined signal for all subjects
combinedSignal = [];
for i = 1:numSubjects
    combinedSignal = [combinedSignal; signal{i}];
end
Y = zeros(1, size(combinedSignal,1));
for i = 1:size(Y,2)
    if any(combinedSignal(i,:))
        Y(i) = 1;
    end
end

%Yhat from group lasso
Yhat = zeros(numLambda, size(combinedSignal,1));
for index = 1:numLambda
    estimatedSignal = [];
    for i = 1:numSubjects
        temp = zeros(size(signal{i}));
        betaIndices = return_index(grpsizes, rois, i);
        idx = [rois{i}{:}];
        for j = 1:numel(idx)
            temp(idx(j),:) = beta{index}(betaIndices(j), :);
        end
        estimatedSignal = [estimatedSignal; temp];
    end
    
    %get Yhat
    for i = 1:size(Yhat,2)
        Yhat(index, i) = norm(estimatedSignal(i,:));
    end
    
    %normalize to [0,1]
    if max(Yhat(index,:) > 0)
        Yhat(index,:) = Yhat(index,:) ./ max(Yhat(index,:));
    end
end

%compute for minnorm solution
idx = [];
for i = 1:numSubjects
    idx = [idx, [rois{i}{:}]];
end
estimatedSignal = zeros(size(estimatedSignal));
for i = 1:numel(idx)
    estimatedSignal(idx(i),:) = betaMinNorm(i,:);
end
YhatMinNorm = zeros(1, size(combinedSignal,1));
for i = 1:numel(YhatMinNorm)
    YhatMinNorm(i) = norm(estimatedSignal(i,:));
end
YhatMinNorm = YhatMinNorm ./ max(YhatMinNorm);
