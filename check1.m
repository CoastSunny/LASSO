clear; clc

load new_data.mat;
rand('state', 0);

alpha = 0;


%% data processing

%input the groups here%%%%%%%%%%%
idx = [3 5];
%Note that in what follows, we index these from 1, i.e. if idx=[7,10], then
%X{1} is the design matrix for group 7. The default is to simulate from
%X{1}.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newROIs.Tess_ndx = newROIs.Tess_ndx(idx);
newROIs.corr = newROIs.corr(idx);
newROIs.name = newROIs.name(idx);


group_ids{1}=[1:100]; group_ids{2}=[101:169];


G = numel(idx); no_groups=G;

n = size(EEG_fwd,1);

p = zeros(1,G);
penalty = p;
X = cell(1,G);
D = X; D_U = X; D_D = X; D_V = X;
signal = X;
betaold = cell(1,G);
for i = 1:G
    p(i) = numel(newROIs.Tess_ndx{i});
    X{i} = EEG_fwd(:,newROIs.Tess_ndx{i});
    a = (eye(n)-1/n*ones(n,n))*X{i};
    penalty(i) = sqrt(trace(newROIs.corr{i}*a'*a));
    signal{i} = zeros(numel(newROIs.Tess_ndx{i}),1);
    D{i} = penalty(i)*chol(inv(newROIs.corr{i}));
    %D{i}=eye(size(D{i}))*p(i);
    betaold{i} = zeros(p(i),1);
end

%simulate signal from X{1}
signal{1} = rand(size(signal{1}));


Y = zeros(n,1);
for i = 1:G
    Y = Y + X{i}*signal{i};
end

%% determine lambda path
lambdaMAX = 0;
for i = 1:G
    lambdaMAX = max(lambdaMAX, norm(D{i}'\X{i}'*Y,2)/(1-alpha));
end

%% fitting
INNER_TOL=1e-6; INNER_ITER_MAX=1;
lambda = lambdaMAX/10;
lam_grp=lambda*(1-alpha); lam_lasso=lambda*alpha;
Maxiter =500;

tk=1; 

w_L = cell(1,G);
w_G = w_L;
for i = 1:G
    [D_U{i},D_D{i},D_V{i}]=svd(full(D{i}));  
    %%% bad matlab function naming. this is not svds
    
    
    w_L{i} = randn(p(i),1); w_L{i}=sign(w_L{i}).*min(abs(w_L{i}),lam_lasso);
    w_G{i} = randn(p(i),1); w_G{i}=lam_grp*w_G{i}/norm(w_G{i},2);
end


Lipsh_vec = zeros(1,G);
for i = 1:G
    s=svd(X{i});
   % Lipsh_vec(i) =normest(X{i})^2;
    Lipsh_vec(i)=s(1)^2;
end

kk=0; betanew = betaold;

obj_vals=zeros(Maxiter*G*INNER_ITER_MAX,1);

time_counter=obj_vals;


xold=zeros(sum(p),1); xnew=xold;

A=[X{1},X{2}];

kk=1;


obj_vals(kk)= 0.5*norm(Y - A*xnew,2)^2;



for iter =1: Maxiter

%iter
obj_val_grp=obj_vals(kk);

for group = 1: no_groups   
     
     xcouple=xold;
     inner_tol=1; inner_iter=0;
     ids=group_ids{group};
     Lipsh_curr=Lipsh_vec(group);
     
  group1=setdiff([1,2],group);   ids1=group_ids{group1};
 
  if (norm(D{group}'\(X{group}'*(Y- X{group1}*xnew(ids1))),2) < lam_grp)
       xnew(ids)=0*xnew(ids); 
     obj_vals(kk) = 0.5*norm(Y - A*xnew,2)^2 + lam_lasso*sum(abs(xnew)) +...
     lam_grp*norm(D{1}*xnew(group_ids{1}),2) + lam_grp*norm(D{2}*xnew(group_ids{2}),2);
    %%xold=xnew;
     xold = xnew;
     xcouple=xnew; 
  
  else      
     
     
     
     while (inner_tol > INNER_TOL)&&(inner_iter<INNER_ITER_MAX)
     inner_iter=inner_iter +1;
     kk=kk+1;
     gradvec= - A'*(Y - A*xcouple);
     gradvec_group=gradvec(ids);
     
     x_0 = xcouple(ids) - gradvec_group/Lipsh_curr; 

%{
     xnew(ids) = prox_map_sparse_grp_lasso(x_0, lam_lasso/Lipsh_curr,lam_grp/Lipsh_curr)  ; 
 obj_vals(kk) = 0.5*norm(Y - A*xnew,2)^2 + lam_lasso*sum(abs(xnew)) +...
     lam_grp*norm(xnew(group_ids{1}),2) + lam_grp*norm(xnew(group_ids{2}),2);
%}

 
        
     t=tic;
[xnew(ids),obj_valsP,w_L{group},w_G{group}]=prox_map_gen_sp_grp_lasso(x_0, D_U{group}, D_D{group},D_V{group} ,lam_lasso/Lipsh_curr,...
    lam_grp/Lipsh_curr, w_L{group},w_G{group},1000) ;

obj_vals(kk) = 0.5*norm(Y - A*xnew,2)^2 + lam_lasso*sum(abs(xnew)) +...
     lam_grp*norm(D{1}*xnew(group_ids{1}),2) + lam_grp*norm(D{2}*xnew(group_ids{2}),2);
time_counter(kk)=toc(t);
%%xold=xnew;
    
     
     inner_tol= norm(xnew - xold)/(norm(xold)+10^-6);
     
     xold = xnew;
     xcouple=xnew;
     
     end

  end
     
end    

obj_vals(kk)

if ( ( - obj_vals(kk) + obj_val_grp)/obj_val_grp < 10^-5)
%    break
end

     xold=xnew;
     

end














