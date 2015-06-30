function [errorCurve, minNormError] = compute_error_curves(Y, beta, betaMinNorm, rois, grpsizes, stackedForwards)

numSubjects = numel(rois);
numLambda = numel(beta);
errorCurve = zeros(1, numLambda);
for i = 1:numLambda
    betahat = [];
    for s = 1:numSubjects
        betahat = [betahat; beta{i}(return_index(grpsizes, rois, s), :)];
    end
    errorCurve(i) = norm(Y-stackedForwards*betahat, 'fro');
end

minNormError = norm(Y-stackedForwards*betaMinNorm, 'fro');