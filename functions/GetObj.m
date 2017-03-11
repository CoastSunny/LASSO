function obj = GetObj(Y, beta, EEG_fwd, ROIs,lam_grp,lam_lasso)

obj = zeros(1,numel(beta));
G = numel(beta{1});
n = numel(Y);
penalty = zeros(1,G);
D = cell(1,G); L = D; X = D;
for i = 1:G
    X{i} = EEG_fwd(:,ROIs.Tess_ndx{i});
    a = (eye(n)-1/n*ones(n,n))*X{i};
    penalty(i) = sqrt(trace(ROIs.corr{i}*a'*a));
    L{i} = 1/penalty(i)*chol(ROIs.corr{i});
    D{i} = L{i}'\speye(size(L{i}));
end


for i = 1:numel(beta)
    res = Y;
    for j = 1:G
        res = res - X{j}*beta{i}{j};
        obj(i) = obj(i) + lam_grp(i)*norm(D{j}*beta{i}{j},2) + lam_lasso(i)*norm(beta{i}{j},1);
    end
    obj(i) = obj(i) + 0.5*norm(res,2)^2;
end