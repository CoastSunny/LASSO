function s = SoftThres(X,r,lam_lasso)

z = X'*r;
s = sign(z).*max(abs(z)-lam_lasso);