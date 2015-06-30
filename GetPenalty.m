
function pen = GetPenalty(D,U,D_sq,sigma_sq,lam_grp,lam_lasso,betaold)

pen = 0;
G_small = numel(D);
G = numel(betaold);
for i = 1:G_small
    pen = pen + lam_grp*norm(D{i}*betaold{i},2) + lam_lasso*norm(betaold{i},1);
end

if G > G_small
    for i = (G_small+1):G
        tmp = U{i}'*betaold{i};
        pen = pen + lam_grp*sqrt(tmp'*(tmp./D_sq{i})+1/sigma_sq(i)*(betaold{i}'*betaold{i}-tmp'*tmp)) + ...
            lam_lasso*norm(betaold{i},1);
    end
end