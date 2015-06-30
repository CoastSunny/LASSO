function [x,obj_valsP,w_L,w_G]=prox_map_gen_sp_grp_lasso_newton(x_0, D,D_U,D_D,D_V, lam_lasso,lam_grp, w_L,w_G,Maxiter) 

%% D= D_U*D_D*D_V';
%% min 0.5*\| x - x_0\|_2^2 + \lam_lasso ||x||_1 + \lam_grp \|Dx\|_2

%% DUAL: min \|D'*w_G + w_L - x_0 \|_2^2 \sbt \|w_G\|_2 \leq \lam_grp , \|w_L\|_\infty \leq \lam_lasso


if (isempty(w_L))
w_L=zeros(n_lasso,1);
w_G=zeros(n_grp,1);
end


dual_tol=1; TOL=1e-5;

if(isempty(Maxiter))
    Maxiter=100;
end
    obj_valsD=zeros(Maxiter,1);
obj_valsP=obj_valsD;
iter=1;

temp = (w_G'*D)';

obj_valsD(iter)= -0.5*norm(w_L + temp - x_0,2)^2 + 0.5*norm(x_0,2)^2 ;

while (dual_tol > TOL)&&(iter<Maxiter) 
% max-threshold w_L
iter=iter+1;
target=x_0 - temp;
w_L=sign(target).*min(abs(target),lam_lasso);

% do a l2-prox map wrt w_G
target=x_0 - w_L;

%% D'= D_V*D_D'*D_U' is the svd of X in \|Y - Xb\|_2 , the required form for fit

 w_G=fit_modified_newton(D_V,D_D',D_U,target,lam_grp,1e-5);

% [w_G,tmp]=fit_modified_nest(D,Dtrans,LipshTrans,target,lam_grp,TOL,w_G);

temp = (w_G'*D)';

obj_valsD(iter)= -0.5*norm(target - temp,2)^2 + 0.5*norm(x_0,2)^2 ;


x= x_0 - temp - w_L; %primal variable

obj_valsP(iter)=0.5*norm(x - x_0,2)^2 + lam_lasso*norm(x,1) + lam_grp*norm(D*x,2);

%%dual_tol=obj_valsP(iter) - obj_valsD(iter);

dual_tol=abs(obj_valsP(iter) - obj_valsP(iter-1))/abs(obj_valsP(iter));


end

obj_valsD=obj_valsD(1:iter);
obj_valsP=obj_valsP(1:iter);




