clear; clc

load Subject_48_initialization

rand('seed', 0);

alpha = 0;
NUM_OF_LAMBDA = 100;

RANK_APPROX = 100;
G_small = 14;

%% data processing

%input the groups here%%%%%%%%%%%
idx = 1:14;
ndx_areas = [3 5];
noise_level = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROIs.Tess_ndx = ROIs.Tess_ndx(idx);
ROIs.corr = ROIs.corr(idx);
% for i = 1:14
%     ROIs.corr{i} = eye(size(ROIs.corr{i}));
% end
ROIs.name = ROIs.name(idx);
ROIs.ndx = ROIs.ndx(idx);
[Y, signal, data, strength] = GenerateData(ROIs, ndx_areas, VertConn, EEG_fwd, noise_level);
Y = Y(:,1);


G = numel(idx);

n = size(EEG_fwd,1);

p = zeros(1,G);
penalty = p;
X = cell(1,G);
D = cell(1,G_small); D_U = D; D_D = D; D_V = D; L = D; U=X; D_sq=X;
sigma_sq = zeros(1,G);
betaold = cell(1,G);
betanew = cell(1,NUM_OF_LAMBDA);
for i = 1:G
    p(i) = numel(ROIs.Tess_ndx{i});
    X{i} = EEG_fwd(:,ROIs.Tess_ndx{i});
    a = (eye(n)-1/n*ones(n,n))*X{i};
    if i <= G_small
        penalty(i) = sqrt(trace(ROIs.corr{i}*a'*a));
        L{i} = 1/penalty(i)*chol(ROIs.corr{i});
        D{i} = L{i}'\speye(size(L{i}));
    else
        [U{i}, D_sq{i}] = eigs(ROIs.corr{i},RANK_APPROX,'LM');
        sigma_sq(i) = trace(ROIs.corr{i}-U{i}*D_sq{i}*U{i}')/(p(i)-RANK_APPROX);
        ROIs.corr{i} = U{i}*D_sq{i}*U{i}'+sigma_sq(i)*(eye(p(i))-U{i}*U{i}');
        penalty(i) = sqrt(trace(ROIs.corr{i}*a'*a));
        D_sq{i} = diag(D_sq{i})/penalty(i)^2;
        sigma_sq(i) = sigma_sq(i)/penalty(i)^2;
    end
    betaold{i} = zeros(p(i),1);
end

for i = 1:NUM_OF_LAMBDA
    betanew{i} = betaold;
end

%% fitting parameters
INNER_TOL=1e-6; INNER_ITER_MAX=1; TOL = 1e-5;
Maxiter =400;
obj_vals = zeros(NUM_OF_LAMBDA,Maxiter*G*INNER_ITER_MAX);
Lipsh_vec = zeros(1,G);
for i = 1:G
    d = svd(X{i}, 'econ');
    Lipsh_vec(i) = d(1)^2;
    if i < 15
        [D_U{i},D_D{i},D_V{i}]=svd(full(D{i}),'econ');
    end
end
%% determine lambda path
lambdaMAX = 0;
for i = 1:G
    if i <= G_small
        lambdaMAX = max(lambdaMAX, norm(L{i}*X{i}'*Y,2)/(1-alpha));
    else
        lambdaMAX = max(lambdaMAX, 1/penalty(i)*sqrt(Y'*X{i}*ROIs.corr{i}*X{i}'*Y)/(1-alpha));
    end
end
lambda = lambdaMAX:(1-lambdaMAX)/(NUM_OF_LAMBDA):1;
lambda = lambda(1:end-1);
%% fitting

lam_grp=lambda*(1-alpha); lam_lasso=lambda*alpha;

mag = zeros(NUM_OF_LAMBDA, G);
dev = zeros(1,numel(lambda));

fprintf('Run\tObjective\n');
fprintf('---------------------------------------\n');
hold
color = zeros(G,3);
for i = 1:G
    color(i,:) = rand(1,3);
end
for runs = 1:numel(lambda)
    
    residual = GetResidual(Y,X,betaold);
    pen = GetPenalty(D,U,D_sq,sigma_sq,lam_grp(runs),lam_lasso(runs),betaold);
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
    while tol > TOL && iter <= Maxiter
        iter = iter + 1;
        obj_old = obj;
        for group = 1:G
            kk = kk + 1;
            partial_res = residual + X{group}*betaold{group};
            cue = 1;
            if group <= G_small
                pen = pen - lam_lasso(runs)*norm(betaold{group},1) -...
                    lam_grp(runs)*norm(D{group}*betaold{group},2);
                if norm(L{group}*SoftThres(X{group},partial_res,lam_lasso(runs)),2) <= lam_grp(runs)
                    cue = 0;
                end
            else
                tmp = U{group}'*betaold{group};
                pen = pen - lam_grp(runs)*sqrt(tmp'*(tmp./D_sq{group})+1/sigma_sq(group)*(betaold{group}'*betaold{group}-tmp'*tmp)) - ...
                    lam_lasso(runs)*norm(betaold{group},1);
                S = SoftThres(X{group},partial_res,lam_lasso(runs));
                tmp = U{group}'*S;
                %if S'*ROIs.corr{group}*S...
                if tmp'*(tmp.*D_sq{group})+sigma_sq(group)*(S'*S-tmp'*tmp) <= lam_grp(runs)^2
                    cue = 0;
                end
            end
            if cue == 0
                betanew{runs}{group} = 0*betanew{runs}{group};
                residual = partial_res;
                obj_vals(runs,kk) = 0.5*norm(residual)^2 + pen;
                betaold{group} = betanew{runs}{group};
            else
                betacouple=betaold{group};
                inner_tol=1; inner_iter=0;
                Lipsh_curr=Lipsh_vec(group);
                while (inner_tol > INNER_TOL)&&(inner_iter<INNER_ITER_MAX)
                    inner_iter=inner_iter +1;
                    gradvec_group= - X{group}'*(partial_res - X{group}*betacouple);
                    beta_0 = betacouple - gradvec_group/Lipsh_curr; 
                    t=tic;
                    if group <= G_small
                        [betanew{runs}{group},obj_valsP,w_L{group},w_G{group}]=prox_map_gen_sp_grp_lasso_newton(beta_0, D{group}, D_U{group},D_D{group},D_V{group}, lam_lasso(runs)/Lipsh_curr,...
                        lam_grp(runs)/Lipsh_curr, w_L{group},w_G{group},400) ;                        
                        pen = pen + lam_lasso(runs)*norm(betanew{runs}{group},1) +...
                            lam_grp(runs)*norm(D{group}*betanew{runs}{group},2);
                    else
                        [betanew{runs}{group},obj_valsP,w_L{group},w_G{group}]=prox_map_gen_sp_grp_lasso_newton_large(beta_0, D_sq{group}, U{group},sigma_sq(group),lam_lasso(runs)/Lipsh_curr,...
                        lam_grp(runs)/Lipsh_curr, w_L{group},w_G{group},400);
                        tmp = U{group}'*betanew{runs}{group};
                        pen = pen + lam_grp(runs)*sqrt(tmp'*(tmp./D_sq{group})+1/sigma_sq(group)*(betanew{runs}{group}'*betanew{runs}{group}-tmp'*tmp)) + ...
                            lam_lasso(runs)*norm(betanew{runs}{group},1);
                    end
                    residual = partial_res - X{group}*betanew{runs}{group};
                    obj_vals(runs,kk) = 0.5*norm(residual)^2 + pen;
                    time_counter(runs,kk)=toc(t);       
                    inner_tol= norm(betanew{runs}{group} - betaold{group})/(norm(betaold{group})+10^-8);
                    betaold{group} = betanew{runs}{group};
                end
            end 
        end     
        obj = obj_vals(runs,kk);
        tol = (obj_old-obj)/obj;    
    end
    fprintf('%d\t%f\n',runs,obj);
    for j = 1:G
        mag(runs,j) = norm(betanew{runs}{j},2)/sqrt(p(j));
    end
    dev(runs) = (Y-residual)'*(Y-residual)/(Y'*Y);
%     a = find(mag(runs,:)>0.2*max(max(mag)));
end

m = find(dev <= 0.8); m = m(end);
[~,m] = max(mag(m,:));
n1 = find(dev <= 0.5); n1 = n1(end);
[~,n1] = max(mag(n1,:));
m = [m n1];
m = unique(m);
m = setdiff(m, ndx_areas);
a = [ndx_areas m];

result = zeros(NUM_OF_LAMBDA, 2);
for i = 1:NUM_OF_LAMBDA
    [result(i,1) result(i,2)] = EvaluateOutput(betanew,i,p,ROIs,signal,VertConn);
end

clf;
% plot(dev, mag(:,a),'LineWidth',2);
YMAX = max(max(mag));
% axis([0 1 0 YMAX]);
Script_create_simulations;
% hold;
% %plot([dev(1) dev(end)],[norm(J_MN) norm(J_MN)],'LineWidth',2,'Color','black');
% legend(ROIs.name(a),'Location','Best');
% text(repmat(dev(end),numel(a),1), mag(end,a), ROIs.name(a), 'FontSize',12);
% xlabel('Variance explained');
% ylabel('Coefficient norm');
% title('Group');grid
% 
% figure;
plot(dev,result,'LineWidth',2);
axis([0 1 0 1]);
hold
plot([dev(1) dev(end)],[AUC AUC],'--','LineWidth',2)
plot([dev(1) dev(end)],[relative_energy relative_energy],'--','LineWidth',2,'Color','green')
xlabel('Variance explained');
legend('AUC','Relative energy','AUC MN','Relative energy MN','Location','Best');
title('Performance measures');
grid

%write Y to file
fid = fopen('Y.txt', 'w');
for i = 1:numel(Y)
    fprintf(fid, '%f\n', Y(i));
end
fclose(fid);