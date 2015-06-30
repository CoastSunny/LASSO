%function xhat =prox_map_gen_sp_grp_lasso(x_0, D_U, D_D,D_V ,lam_lasso,lam_grp, w_L,w_G ) 

%% D= D_U*D_D*D_V';
%% min 0.5*\| x - x_0\|_2^2 + \lam_lasso ||x||_1 + \lam_grp \|Dx\|_2

%% DUAL: min \|D'*w_G + w_L - x_0 \|_2^2 \sbt \|w_G\|_2 \leq \lam_grp , \|w_L\|_\infty \leq \lam_lasso


%% generate data:
clear; clc;
mm=10; nn=50; lam_lasso=1; lam_grp=1;
D=randn(mm,nn); [D_U,D_D,D_V]=svd(D);
w_L=randn(nn,1); w_L=sign(w_L).*min(abs(w_L),lam_lasso);
w_G=randn(mm,1); w_G=lam_grp*w_G/norm(w_G,2);

x_0= randn(nn,1) + D'*w_G + w_L ;
lam_lasso=.5; lam_grp=.5;

n_lasso=length(x_0); 
n_grp=size(D_U,1);

D=D_U*D_D*D_V';

if (isempty(w_L))
w_L=zeros(n_lasso,1);
w_G=zeros(n_grp,1);
end


dual_tol=1; TOL=10^-6;

Maxiter=100;
obj_valsD=zeros(Maxiter,1);
obj_valsP=obj_valsD;
iter=0;

while (dual_tol > TOL)&&(iter<Maxiter) 
% max-threshold w_L
iter=iter+1;
target=x_0 - D'*w_G;
w_L=sign(target).*min(abs(target),lam_lasso);

% do a l2-prox map wrt w_G
target=x_0 - w_L;

%% D'= D_V*D_D'*D_U' is the svd of X in \|Y - Xb\|_2 , the required form for fit

w_G=fit_modified(D_V,D_D',D_U,target,lam_grp,TOL);

obj_valsD(iter)= -0.5*norm(target - D'*w_G,2)^2 +0.5*norm(x_0,2)^2 ;

x= x_0 - D'*w_G - w_L; %primal variable

obj_valsP(iter)=0.5*norm(x - x_0,2)^2 + lam_lasso*norm(x,1) + lam_grp*norm(D*x,2);

dual_tol=obj_valsP(iter) - obj_valsD(iter);


end