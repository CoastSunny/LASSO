clear; clc

load Subject_48_initialization

rand('state', 0);

alpha = 0;
NUM_OF_LAMBDA = 100;

%% data processing

%input the groups here%%%%%%%%%%%
idx = 1:14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROIs.Tess_ndx = ROIs.Tess_ndx(idx);
ROIs.corr = ROIs.corr(idx);
ROIs.name = ROIs.name(idx);
ROIs.ndx = ROIs.ndx(idx);

G = numel(idx);

n = size(EEG_fwd,1);

p = zeros(1,G);
penalty = p;
X = cell(1,G);
D = X; Dtrans = X; L = X;
betaold = cell(1,G);
betanew = cell(1,NUM_OF_LAMBDA);
for i = 1:G
    p(i) = numel(ROIs.Tess_ndx{i});
    X{i} = EEG_fwd(:,ROIs.Tess_ndx{i});
    a = (eye(n)-1/n*ones(n,n))*X{i};
    penalty(i) = sqrt(trace(ROIs.corr{i}*a'*a));
	L{i} = 1/penalty(i)*chol(ROIs.corr{i});
	D{i} = L{i}'\speye(size(L{i}));
    Dtrans{i} = D{i}';
    betaold{i} = zeros(p(i),1);
end

for i = 1:NUM_OF_LAMBDA
    betanew{i} = betaold;
end

[Y, Source_distribution] = GenerateData(ROIs, VertConn, EEG_fwd);
Y = Y(:,1);

%% fitting parameters
INNER_TOL=1e-6; INNER_ITER_MAX=1; TOL = 1e-5;
Maxiter =4000;
obj_vals = zeros(NUM_OF_LAMBDA,Maxiter*G*INNER_ITER_MAX);
Lipsh_vec = zeros(1,G);
LipshTrans = zeros(1,G);
for i = 1:G
    Lipsh_vec(i) = normest(X{i})^2;
    LipshTrans(i) = my_normest(Dtrans{i},D{i})^2;
end
%% determine lambda path
lambdaMAX = 0;
for i = 1:G
    lambdaMAX = max(lambdaMAX, norm(L{i}*X{i}'*Y,2)/(1-alpha));
end
lambda = lambdaMAX:(1-lambdaMAX)/(NUM_OF_LAMBDA-1):1;
%% fitting

lam_grp=lambda*(1-alpha); lam_lasso=lambda*alpha;

mag = zeros(NUM_OF_LAMBDA, G);

fprintf('Run\tObjective\n');
fprintf('---------------------------------------\n');

for runs = 1:numel(lambda)
    
    residual = GetResidual(Y,X,betaold);
    pen = GetPenalty(D,lam_grp(runs),lam_lasso(runs),betaold);
    obj = 0.5*norm(residual)^2 + pen;
    
    tk=1; 

    w_L = cell(1,G);
    w_G = w_L;
    for i = 1:G
        w_L{i} = randn(p(i),1); w_L{i}=sign(w_L{i}).*min(abs(w_L{i}),lam_lasso(runs));
        w_G{i} = randn(p(i),1); w_G{i}=lam_grp(runs)*w_G{i}/norm(w_G{i},2);
    end

    kk=0;
    
    time_counter=obj_vals(runs,:);
    iter = 0;
    tol = 1;
    while tol > TOL 
        iter = iter + 1;
        obj_old = obj;
        for group = 1:G
            kk = kk + 1;
            partial_res = residual + X{group}*betaold{group};
            pen = pen - lam_lasso(runs)*norm(betaold{group},1) -...
            lam_grp(runs)*norm(D{group}*betaold{group},2);
            if norm(L{group}*SoftThres(X{group},partial_res,lam_lasso(runs)),2) <= lam_grp(runs)
                betanew{runs}{group} = 0*betanew{runs}{group};
                residual = partial_res;
                obj_vals(runs,kk) = 0.5*norm(residual)^2 + pen;
            else
                betacouple=betaold{group};
                inner_tol=1; inner_iter=0;
                Lipsh_curr=Lipsh_vec(group);
                while (inner_tol > INNER_TOL)&&(inner_iter<INNER_ITER_MAX)
                    inner_iter=inner_iter +1;
                    gradvec_group= - X{group}'*(partial_res - X{group}*betacouple);
                    beta_0 = betacouple - gradvec_group/Lipsh_curr; 
                    t=tic;
                    [betanew{runs}{group},obj_valsP,w_L{group},w_G{group}]=prox_map_gen_sp_grp_lasso(beta_0, D{group}, Dtrans{group}, LipshTrans(group), lam_lasso(runs)/Lipsh_curr,...
                    lam_grp(runs)/Lipsh_curr, w_L{group},w_G{group},400) ;
                    residual = partial_res - X{group}*betanew{runs}{group};
                    pen = pen + lam_lasso(runs)*norm(betanew{runs}{group},1) +...
                    lam_grp(runs)*norm(D{group}*betanew{runs}{group},2);
                    obj_vals(runs,kk) = 0.5*norm(residual)^2 + pen;
                    time_counter(runs,kk)=toc(t);       
                    inner_tol= norm(betanew{runs}{group} - betaold{group})/(norm(betaold{group})+10^-6);
                    betaold{group} = betanew{runs}{group};
                end
            end
        end     
        obj = obj_vals(runs,kk);
        tol = (obj_old-obj)/obj;     
    end
    fprintf('%d\t%f\n',runs,obj);
end

for i = 1:NUM_OF_LAMBDA
    for j = 1:G
        mag(i,j) = norm(betanew{i}{j},2);
    end
end

plot(mag);
legend(ROIs.name,'Location','NorthWest');


AUC = zeros(1,100);
AUC_close = AUC; AUC_far = AUC; n_MSE = AUC; n_DF = AUC; relative_energy = AUC;
for i = 1:100
    beta = [];
    for j = 1:G
     beta = [beta;betanew{i}{j}];
    end
    beta = [beta; zeros(size(EEG_fwd,2)-numel(beta),1)];
    [ AUC(i) , AUC_close(i) , AUC_far(i) , n_MSE(i) , n_DF(i) , relative_energy(i) ] = Inverse_perf_estimation( Source_distribution , beta , VertConn );
end
