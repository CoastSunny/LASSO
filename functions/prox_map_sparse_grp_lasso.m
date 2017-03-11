function xhat =prox_map_sparse_grp_lasso(x_0, lam_lasso,lam_grp) 

%% min 0.5*\| x - x_0\|_2^2 + \lam_lasso ||x||_1 + \lam_grp \|x\|_2

% soft-threshold x_0

x_0new=sign(x_0).*max(abs(x_0) -lam_lasso,0);

% do a l2-prox map on x_0new

xhat = x_0new*max(1 - lam_grp/norm(x_0new,2),0);


end