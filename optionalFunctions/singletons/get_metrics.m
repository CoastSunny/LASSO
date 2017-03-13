function [auc, relativeEnergy] = get_metrics(beta, source, VertConn, indices)

T = size(beta,2);

auc = zeros(1,T);
relativeEnergy = zeros(1,T);

idx = [indices{:}];
estimatedSource = zeros(size(source));
for t = 1:T
    for i = 1:numel(idx)
        estimatedSource(idx(i),t) = beta(i,t);
    end
    [auc(t),~,~,~,~,relativeEnergy(t)] = Inverse_perf_estimation(source(:,t),estimatedSource(:,t),VertConn);
end

auc = mean(auc);
relativeEnergy = mean(relativeEnergy);