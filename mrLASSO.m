clear all;
close all;

% ADD NECESSARY FUNCTIONS
addpath(genpath(pwd));

% MAIN PARAMETERS
numSubs = 10; % if left empty, use all subs
doMinNorm = true;    % compute minimum norm or load in previous version?
newSubs = false;      % subjects and ROIs used in original paper or new ones
simulateData = true; % simulated or real data?
%projectDir = '~/Dropbox/ONGOING/LASSO/forward_data'; %'/Volumes/svndl/4D2/kohler/SYM_16GR/SOURCE';
projectDir = '/Users/ales/Downloads/forward_data';
condNmbr = 15; % which condition to plot, only relevant if new subs and real data

if ~newSubs;
    if ~simulateData
        error('You cannot run old subjects with real data');
    else
    end
    condNmbr = [];
else
end


% ADDITIONAL PARAMETERS?
stackedForwards = [];
allSubjForwards = {};
ridgeSizes = zeros(1, numSubs);
numComponents = 5;
numCols = 2;
nLambdaRidge = 50;
nLambda = 30; % always end picking the max? (~120)
alphaVal = 1.0817e4;
MAX_ITER = 1e6;

%% GET ROIs AND EXCLUDE MISSING
if newSubs
    tic
    subjectList = subfolders([projectDir,'/*']);
    if isempty(numSubs)
        numSubs = length(subjectList);
    else
    end
    for s=1:numSubs
        ROIs{s} = sl_roiMtx(subjectList{s});
        if s==1
            numROIs = size(ROIs{s}.ndx,2);
            roiNames = ROIs{s}.name;
        else
        end
        roiIdx{s} = ROIs{s}.ndx;
        nanROIs(s,:) = cell2mat(arrayfun(@(x) ~isempty(roiIdx{s}{x}),1:numROIs,'uni',0));
    end
    numROIs = length(find(sum(nanROIs,1)>=numSubs));
    roiIdx = arrayfun(@(x) roiIdx{x}(sum(nanROIs,1)>=numSubs),1:numSubs,'uni',0); % modify roiIdx
    
    for s=1:numSubs
        ROIs{s}.ndx = ROIs{s}.ndx(sum(nanROIs,1)>=numSubs); % modify ROIs
        ROIs{s}.corr = ROIs{s}.corr(sum(nanROIs,1)>=numSubs);
        ROIs{s}.name = ROIs{s}.name(sum(nanROIs,1)>=numSubs);
    end
    
    %% MAKE FORWARD MATRIX, XLIST AND VLIST
    for s=1:numSubs
        fwdMatrix = sl_getFwd(projectDir,subjectList{s});
        ridgeSizes(s) = numel(cat(2,roiIdx{s}{:})); % total ROI size by subject
        stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[roiIdx{s}{:}]));
        
        % make Xlist and Vlist
        % get the first x principle components of the part of the fwdMatrix corresponding to each ROI
        xList{s} = cell(1,numROIs); vList{s} = cell(1,numROIs);
        [xList{s}(:),vList{s}(:)] = arrayfun(@(x) get_principal_components(fwdMatrix(:, roiIdx{s}{x}),numComponents),1:numROIs,'uni',0);
    end
    toc
else
    subjectList = [1,3,4,9,17,35,36,37,39,44,48,50,51,52,53,54,55,66,69,71,75,76,78,79,81];
    if isempty(numSubs)
        numSubs = length(subjectList);
    else
    end
    for s=1:numSubs
        structure = load([projectDir,'/ROI_correlation_many_subjects/ROI_correlation_subj_' num2str(subjectList(s))]);
        ROIs{s} = structure.ROIs;
        roiIdx{s} = ROIs{s}.ndx;
        numROIs = length(roiIdx{1}); % all subjects should have same number of ROIs
        structure = load([projectDir,'/forwardAndRois-skeri' num2str(subjectList(s))]);
        fwdMatrix = structure.fwdMatrix;
        ridgeSizes(s) = numel(cat(2,roiIdx{s}{:})); % total ROI size by subject
        stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[roiIdx{s}{:}]));
        allSubjForwards{s} = fwdMatrix;
        % make Xlist and Vlist
        % get the first x principle components of the part of the fwdMatrix corresponding to each ROI
        xList{s} = cell(1,numROIs); vList{s} = cell(1,numROIs);
        [xList{s}(:),vList{s}(:)] = arrayfun(@(x) get_principal_components(fwdMatrix(:, roiIdx{s}{x}),numComponents),1:numROIs,'uni',0);
    end
end

%% COMPUTE INVERSES
% get number of subjects and conditions
% (should be same for all subjects)
for s=1:numSubs
    % get data
    if simulateData
        if s == 1 % only randomize phase for first subject
            phase = randi([2,10], 1);
        else
        end
        if newSubs
            signalROIs = {'func_V2v-L','func_V4-R'}; % simulate signal coming from these ROIs
        else
            signalROIs = {'V2v-L','V4-R'};
        end
        idx = cell2mat(arrayfun(@(x) cellfind(lower(ROIs{s}.name),lower(signalROIs{x})),1:length(signalROIs),'uni',false));
        if numel(idx)<numel(signalROIs) % check if ROIs exist
            error('Signal ROIs do not exist!')
        else
        end
        SNR = 0.1;
        initStrct = load([projectDir,'/Subject_48_initialization']);
        [Y(128*(s-1)+1:128*s,:), source{s}, signal{s}, noise] = GenerateData(ROIs{s},idx,initStrct.VertConn,allSubjForwards{s},SNR, phase);
    else
        subjId = subjectList{s};
        exportFolder = subfolders(fullfile(projectDir,subjId,'Exp_MATL_*'),1);
        exportFileList = subfiles(fullfile(exportFolder{1},'Axx_c*.mat'));
        Axx = load([exportFolder{1},'/',exportFileList{condNmbr}]);
        Y(size(Axx.Wave,2)*(s-1)+1:size(Axx.Wave,2)*s,:) = Axx.Wave';
    end
end

%  GENERATE X AND V
X = []; V = [];
for g = 1:numROIs
    tempX = [];
    tempV = [];
    for s = 1:numSubs
        tempX = blkdiag(tempX, xList{s}{g});
        tempV = blkdiag(tempV, vList{s}{g}(:,1:numComponents));
    end
    X = [X, tempX];
    V = blkdiag(V, tempV);
end
grpSizes = numComponents*numSubs*ones(1,numROIs);
indices = get_indices(grpSizes);
penalties = get_group_penalties(X, indices);

% center Y, X, stackedForwards
%scal isn't a matlab builtin, changed to use bsxfun() --jma
%Y = scal(Y, mean(Y));
%X = scal(X, mean(X));

Y = bsxfun(@minus,Y, mean(Y));
X = bsxfun(@minus,X, mean(X));

%changed to bsxfun --jma
%stackedForwards = scal(stackedForwards, mean(stackedForwards));
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));

n = numel(Y);
ssTotal = norm(Y, 'fro')^2 / n;

% minumum norm solution
if doMinNorm
    disp('Generating minimum norm solution');
    [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Y, nLambdaRidge);
    
    
    %This code does the min-norm how mrCurrent does it, So you can compare
    %betaComp with betaMinNorm to see if the results are comperable.
    %this code is basically what's in minimm_norm() but wanted a
    %doublecheck
    [u,s,v] = csvd(stackedForwards);
    lambda = gcv(u,s,Y,'Tikh',100);  
    %Tikhonov regularized inverse matrix
    reg_s = diag( s ./ (s.^2 + lambda^2 ));
    sol = v * reg_s * u';
    betaComp = sol*Y;
    lambdaMinNorm = lambdaMinNorm^2;
    rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotal;
    if simulateData
        save('~/Desktop/sim_minNorm.mat','betaMinNorm','lambdaMinNorm','rsquaredMinNorm');
    else
        save('~/Desktop/minNorm.mat','betaMinNorm','lambdaMinNorm','rsquaredMinNorm');
    end
else
    disp('Loading minimum norm solution');
    if simulateData
        load('~/Desktop/sim_minNorm.mat');
    else
        load('~/Desktop/minNorm.mat');
    end
end

% use first 2 columns of v as time basis
% PK: moved this down to avoid conflict with v generated above
[~, ~, v] = svd(Y);
Ytrans = Y * v(:, 1:numCols);
%Ytrans = scal(Ytrans, mean(Ytrans));
Ytrans = bsxfun(@minus,Ytrans, mean(Ytrans));

% sequence of lambda values
lambdaMax = max(cell2mat(arrayfun(@(x) norm(X(:,indices{x})'*Ytrans, 'fro')/penalties(x),1:numROIs,'uni',0)));
lambdaMax = lambdaMax + 1e-4;
lambdaGrid = lambdaMax * (0.01.^(0:1/(nLambda-1):1));
tol = min(penalties) * lambdaGrid(end) * 1e-5;
if alphaVal > 0
    tol = min([tol, 2*alphaVal*1e-5]);
end

ridgeRange = [0 cumsum(ridgeSizes)];
roiSizes = zeros(1,numROIs); %total size of each region summed over all subjects
for g = 1:numROIs
    roiSizes(g) = sum(cell2mat(arrayfun(@(x) numel(roiIdx{x}{g}),1:numSubs,'uni',0)));
end

% OLS FIT
betaOls = (X'*X + alphaVal*eye(size(X,2))) \ (X'*Ytrans);

% FITTING
betaInit = zeros(size(X,2), numCols);
betaVal = cell(1, nLambda);
objValues = cell(1, nLambda);
gcvError = zeros(1, nLambda);
df = zeros(1, nLambda);
indexer = cell2mat(arrayfun(@(x) return_index(roiSizes, roiIdx, x), 1:numSubs,'uni',0));
for i = 1:nLambda
    [betaVal{i}, objValues{i}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(i), alphaVal, tol, MAX_ITER, penalties, indices);
    betaInit = betaVal{i};
    betaVal{i} = V * betaVal{i} * v(:,1:numCols)'; %transform back to original space (permuted forward matrices)
    rss = norm(Y-stackedForwards*betaVal{i}(indexer, :), 'fro')^2 / n;
    [gcvError(i), df(i)] = compute_gcv(rss, betaInit, betaOls, grpSizes, n);
end
[~, bestIndex] = min(gcvError);

% compute average of average metrics
for s = 1:numSubs
    range = cell2mat(arrayfun(@(x)  numel(roiIdx{s}{x}),1:numROIs,'uni',false));
    range = [0 cumsum(range)];
    tempMinNorm = betaMinNorm(ridgeRange(s)+1:ridgeRange(s+1), :);
    temp = betaVal{bestIndex}(return_index(roiSizes, roiIdx, s), :);
    regionActivityMinNorm(:,:,s) = cell2mat(arrayfun(@(x) mean(tempMinNorm(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionActivity(:,:,s) = cell2mat(arrayfun(@(x) mean(temp(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
end
regionActivityMinNorm = mean(regionActivityMinNorm,3);
regionActivity = mean(regionActivity,3);


%% MAKE FIGURE
    
close all
leftIdx = cell2mat(arrayfun(@(x) ~isempty(strfind(ROIs{1}.name{x},'-L')), 1:length(ROIs{1}.name),'uni',false));
fontSize = 12;
gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Arial','linewidth',1};

tempColors = colormap(jet);
tempIdx = round(linspace(0,length(tempColors),length(find(leftIdx==1))+1));
roiColors = tempColors(tempIdx(2:end),:);



figH(1) = figure;
subplot(2,1,1);
hold on;
selectedRois = find(leftIdx==1);

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(regionActivityMinNorm(selectedRois,:)','linewidth',2)

xlim([0,size(regionActivityMinNorm,2)]);
legend(ROIs{1}.name(leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;

subplot(2,1,2);
hold on;
selectedRois = find(leftIdx==0); %Select left ROI responses

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(regionActivityMinNorm(selectedRois,:)','linewidth',2)


xlim([0,size(regionActivityMinNorm,2)]);
legend(ROIs{1}.name(~leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;
title('minNorm')

figH(2) = figure;
subplot(2,1,1);
hold on;
selectedRois = find(leftIdx==1);

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(regionActivity(selectedRois,:)','linewidth',2)

xlim([0,size(regionActivity,2)]);
legend(ROIs{1}.name(leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;
subplot(2,1,2);
hold on;
selectedRois = find(leftIdx==0);

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(regionActivity(selectedRois,:)','linewidth',2)

legend(ROIs{1}.name(~leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off; 
pos = get(figH(1), 'Position');
pos(3) = pos(3)*2; % Select the height of the figure in [cm]
set(figH, 'Position', pos);
title('lasso');
