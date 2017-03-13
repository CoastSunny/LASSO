%pads beta with zeros, and then calls Inverse_perf_estimation. This should be fair since both the minimum norm and group lasso solutions are obtained only from the visual areas.

function [auc, aucClose, aucFar, mse, n_DF, relativeEnergy] = get_metrics(beta, source, VertConn, indices)

T = size(beta,2);

auc = zeros(1,T);
aucClose = zeros(1,T);
aucFar = zeros(1,T);
mse = zeros(1,T);
n_DF = zeros(1,T);
relativeEnergy = zeros(1,T);

idx = [indices{:}];
estimatedSource = zeros(size(source));
for t = 1:T
    for i = 1:numel(idx)
        estimatedSource(idx(i),t) = beta(i,t);
    end
    [auc(t),aucClose(t), aucFar(t), mse(t), n_DF(t),relativeEnergy(t)] = Inverse_perf_estimation(source(:,t),estimatedSource(:,t),VertConn);
end

auc = mean(auc);
aucClose = mean(aucClose);
aucFar = mean(aucFar);
mse = mean(mse);
n_DF = mean(n_DF);
relativeEnergy = mean(relativeEnergy);
