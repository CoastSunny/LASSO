function [AUC, relative_energy] = EvaluateOutput(betanew, path, p, ROIs, signal, VertConn)

beta = zeros(sum(p), 100);
a = [0 cumsum(p)];
for i = 1:100
    for j = 1:14
        beta(a(j)+1:a(j+1),i) = betanew{i}{j};
    end
end

idx = [ROIs.ndx{1:numel(p)}];

[ AUC , ~ , ~ , ~ , ~ , relative_energy ] = Inverse_perf_estimation( signal(idx) , beta(:,path) , VertConn );