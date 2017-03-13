%%
clear all;
close all;
disp('mrLASSO - Beginning Data Analysis')
disp('Setting up directory structure and some parameters')


msgTxt = {'This script implements fMRI ROI informed group constrained source localization',...
    '1) Compute minimum-norm inverse solution',...
    '2) Compute the Group LASSO inverse solution',...
    '3) Plot waveforms from ROI',...
    '4) Plot topographies',...
    'Press OK to continue'};

uiwait(msgbox(msgTxt,'Information','modal'));
 

% Construct a questdlg with three options
choice = questdlg(['Would you like to Recalculate the minimum norm inverse solution?'...
    ' This calculuation may take 20 minutes'], ...
    'Calculate MinNorm Inverse?',...
    'Recalculate Inverse Solution', ...
    'Use Precalculated Inverse Solution','Use Precalculated Inverse Solution');
switch choice
    case 'Recalculate Inverse Solution', 
        doMinNorm = true;
    case 'Use Precalculated Inverse Solution',
        doMinNorm = false;
end

% Setup paths and directories
myLocation = mfilename('fullpath');
myDir = fileparts(myLocation);
addpath(genpath(myDir));

datafileDir = fullfile(myDir,'datafiles');
if ~exist(datafileDir,'dir')
    uiwait(msgbox('Could not find datafile directory, please select it','Information','modal'));
    datafileDir = uigetdir(myDir,'Pick directory containing data')
end


setpref('mrLASSO','scalpFileDir',fullfile(datafileDir,'anatomy'));


% MAIN PARAMETERS
numSubs = 10; % if left empty, use all subs
simulateData = false; % simulated or real data?
preselectSubjects = false; 
%projectDir = '~/Dropbox/ONGOING/LASSO/forward_data'; %'/Volumes/svndl/4D2/kohler/SYM_16GR/SOURCE';
%projectDir = '/Users/ales/Downloads/forward_data';
projectDir = fullfile(datafileDir,'eegdata');
forwardDir = fullfile(datafileDir,'forwardSolutions');
condNmbr = 8; % which condition to plot,

% ADDITIONAL PARAMETERS?
stackedForwards = [];
allSubjForwards = {};
ridgeSizes = zeros(1, numSubs);
numComponents = 3;
numCols = 5;
lowPassNF1 = true;

nLambdaRidge = 50;
nLambda = 30; % always end picking the max? (~120)
alphaVal = 1.0817e4;
MAX_ITER = 1e6;

%% GET ROIs AND EXCLUDE PARTICIPANTS WITH MISSING INFO
disp('Loading participant ROI information.');

dirList = dir(projectDir);

dirNames = {dirList([dirList(:).isdir]).name};
dirNames = {dirNames{~strncmp('.',dirNames,1)}};
numSubs = length(dirNames);

for iSubj=1:numSubs,
    subjectList(iSubj) = str2num(dirNames{iSubj}(end-3:end));
end

%     numSubs = ;
%     subjectList = subjectList(1:numSubs);
%     dirNames = dirNames(1:numSubs);
    
    if isempty(numSubs)
        numSubs = length(subjectList);
    else
    end
    for s=1:numSubs
        roiFile = fullfile(forwardDir,'ROI_correlation_many_subjects',...
            ['ROI_correlation_subj_' num2str(subjectList(s))]);
        structure = load(roiFile);
        ROIs{s} = structure.ROIs;
        roiIdx{s} = ROIs{s}.ndx;
        numROIs = length(roiIdx{1}); % all subjects should have same number of ROIs
        fwdFile = fullfile(forwardDir,['forwardAndRois-skeri' num2str(subjectList(s))]);
        structure = load(fwdFile);
        fwdMatrix = structure.fwdMatrix;
        ridgeSizes(s) = numel(cat(2,roiIdx{s}{:})); % total ROI size by subject
        stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[roiIdx{s}{:}]));
        allSubjForwards{s} = fwdMatrix;
        % make Xlist and Vlist
        % get the first x principle components of the part of the fwdMatrix corresponding to each ROI
        xList{s} = cell(1,numROIs); vList{s} = cell(1,numROIs);
        [xList{s}(:),vList{s}(:)] = arrayfun(@(x) get_principal_components(fwdMatrix(:, roiIdx{s}{x}),numComponents),1:numROIs,'uni',0);
    end


%% Load the EEG data


disp('Loading EEG data')

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
        initFile = fullfile(forwardDir,'Subject_48_initialization');
        initStrct = load(initFile);
        
        [Y(128*(s-1)+1:128*s,:), source{s}, signal{s}, noise] = GenerateData(ROIs{s},idx,initStrct.VertConn,allSubjForwards{s},SNR, phase);
    else
        subjId = dirNames{s};
        exportFolder = subfolders(fullfile(projectDir,subjId,'Exp_MATL_*'),1);
        exportFileList = subfiles(fullfile(exportFolder{1},'Axx_c*.mat'));
        Axx = load([exportFolder{1},'/',exportFileList{condNmbr}]);
        
        if lowPassNF1
            num2keep = 5;
            nf1 = Axx.i1F1;
            axxIdx = (nf1:nf1:10*nf1)+1;
            dftIdx = (1:10)+1;
            dft = dftmtx(size(Axx.Wave,1));
            sinWaveforms = imag(dft);
            cosWaveforms = real(dft);
            wave = cosWaveforms(:,dftIdx)*Axx.Cos(axxIdx,:)-sinWaveforms(:,dftIdx)*Axx.Sin(axxIdx,:);
            Axx.Wave = wave;
        end
        
    
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
Y = bsxfun(@minus,Y, mean(Y));
X = bsxfun(@minus,X, mean(X));
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));

n = numel(Y);
[u1, s1, v1] = svd(Y);
Ylo = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
ssTotal = norm(Ylo, 'fro')^2 / n;

%% Find Minimum Norm Solution
% minumum norm solution
if doMinNorm
    disp('Generating minimum norm solution');
    [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Ylo, nLambdaRidge);
    
    
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
    YhatMN=stackedForwards*betaMinNorm;
    
    saveDir = datafileDir;
    
    if simulateData
        save(fullfile(saveDir,'simulation_precalculatedMinNorm.mat'),'betaMinNorm','lambdaMinNorm','rsquaredMinNorm','YhatMN');
    else
        save(fullfile(saveDir,'precalculatedMinNorm.mat'),'betaMinNorm','lambdaMinNorm','rsquaredMinNorm','YhatMN');
    end
else
    disp('Loading minimum norm solution');
      saveDir = datafileDir;
    
    if simulateData
        load(fullfile(saveDir,'simulation_precalculatedMinNorm.mat'));
    else
        load(fullfile(saveDir,'precalculatedMinNorm.mat'));
    end
end


%% Group LASSO
% use first 2 columns of v as time basis
% PK: moved this down to avoid conflict with v generated above

disp('Starting group LASSO');
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
YhatLASSO = stackedForwards*betaVal{bestIndex}(indexer, :);

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

baselineIdx = 700:780;
regionActivity = bsxfun(@minus,regionActivity,mean(regionActivity(:,baselineIdx),2));
regionActivityMinNorm = bsxfun(@minus,regionActivityMinNorm,mean(regionActivityMinNorm(:,baselineIdx),2));



%% MAKE FIGURES
    
msgTxt = 'Inverse solutions completed. Now making plots';

disp(msgTxt)



%% Make Figure 7 butterfly plot
msgTxt = 'Recreating Figure 7 butterfly plot';
disp(msgTxt)

figure;
t = 0:Axx.dTms:(Axx.nT-1)*Axx.dTms;
unstackedData = reshape(Y,128,numSubs,780);
grandMeanData = squeeze(mean(unstackedData,2));

plot(t,grandMeanData,'k')
xlabel('Time (ms)')
ylabel('Potential (microvolts)')

%% Make Figure 9: 

msgTxt = 'Recreating Figure 9: Plotting waveforms from ROIs in separate subplots';
disp(msgTxt)

scale = 1e3;
fontSize = 12;
gcaOpts = {'tickdir','out','box','off','fontsize',22,'fontname','Arial',...
    'linewidth',3,'YLim',scale*[-2 2]*1e-4,'clipping','off','XLim',[0 1000],...
    'ticklength',[0.0300 0.050]};
leftIdx = cell2mat(arrayfun(@(x) ~isempty(strfind(ROIs{1}.name{x},'-L')), 1:length(ROIs{1}.name),'uni',false));

selectedRois = find(leftIdx==1);
nRoi=length(selectedRois);
%% Min Norm
figH(1) = figure;


for iRoi = 1:nRoi,
    subplot(ceil(nRoi/2),2,iRoi);
%    subplot(nRoi,1,iRoi);
    hold on;
    
    %For plotting left hemisphere activation
    selectedRois = find(leftIdx==1);
    %hemiName = 'left';
    plot(t,scale*regionActivityMinNorm(selectedRois(iRoi),:)','k-','linewidth',3)

    %For plotting right hemisphere activation
    selectedRois = find(leftIdx==0);
    % hemiName = 'right';
    plot(t,scale*regionActivityMinNorm(selectedRois(iRoi),:)','k:','linewidth',3)
    
    xlim([0,size(regionActivityMinNorm,2)]);    
    set(gca,gcaOpts{:})
    if iRoi == 1;
        title({'Minimum Norm';'Separated waveforms'})
    end
   
    if iRoi ~= nRoi,
        set(gca,'xticklabel','','yticklabel','');
        
    end
    thisRoiName =ROIs{1}.name{selectedRois(iRoi)};
    thisRoiName = thisRoiName(1:end-2);
    text(-320, 0,thisRoiName,'fontsize',22)
    
    axis on
    hline = refline(0,0);
    hline.Color = 'k';
end


%Set the figure position
pos = [0 0 600 900];
set(figH, 'Position', pos);
legend('Left','Right')


%% Lasso
figH(1) = figure;
hold on;


for iRoi = 1:nRoi,
    subplot(ceil(nRoi/2),2,iRoi);%  subplot(nRoi,1,iRoi);
    hold on;
    
    %For plotting left hemisphere activation
    selectedRois = find(leftIdx==1);
    % hemiName = 'Left';
    plot(t,scale*regionActivity(selectedRois(iRoi),:)','k-','linewidth',3)
    
    %For plotting right hemisphere activation
    selectedRois = find(leftIdx==0);
    % hemiName = 'right';
    plot(t,scale*regionActivity(selectedRois(iRoi),:)','k:','linewidth',3)
    
    xlim([0,size(regionActivity,2)]);    
    set(gca,gcaOpts{:})
    if iRoi == 1;
        
        title({'Group Lasso';'Separated waveforms'})
    end
   
    if iRoi ~= nRoi,
        set(gca,'xticklabel','','yticklabel','');
        
    end
    thisRoiName =ROIs{1}.name{selectedRois(iRoi)};
    thisRoiName = thisRoiName(1:end-2);
    text(-320, 0,thisRoiName,'fontsize',22)
    
    axis on
    hline = refline(0,0);
    hline.Color = 'k';
end

%Set the figure position
pos = [600 0 600 900];
set(figH, 'Position', pos);

legend('Left','Right')


%% Topo Reconstructions

%Because the topographies are overlayed on high-resolution scalps
%in order to plot lots of them without performance degradation we plot one
%at a time and take a snapshot of the figure before destroying the figure
%and rendering a new scalp.  For 3D view render a single participant

msgtxt = [ 'Recreating Figure 10'];
disp(msgtxt);

% Construct a questdlg with three options
plotChoice = questdlg(['Would you like to plot topographies?'...
    'The figures are somewhat slow to render'], ...
    'Plot Topographies?',...
    'Yes', ...
    'No','Yes');
switch plotChoice
    case 'Yes',
        
        numSubs = size(unstackedData,2);
        
        %Reduced dimension data (Ylo) that we fit.
        unstackedData = reshape(Ylo,128,numSubs,780);
        %MinNorm Solution
        unstackedYhatMN = reshape(YhatMN,128,numSubs,780);
        %LASSO Solution
        unstackedYhatLASSO = reshape(YhatLASSO,128,numSubs,780);

        numChoice = questdlg({'How many participants to plot?';...
            '1) Single is a 3D rotatable render'; ...
            '2) All is a snapshot of all scalps on one figure'},...
            'Plot Topographies?',...
            'Single', ...
            'All','Single');
        switch numChoice
            case 'Single'
                subjList = 6;
                closeTopoFigs = false;
            case 'All'
                subjList = 1:numSubs;
                closeTopoFigs = true;
        end
            
            %Time index to plot. 196th sample = 251.3 ms.
iT = 196;

%Set a default position for the plots
pos =  [80   300   362   490];


%Create Matrices to save views of topographies.
%Basically we take snapshots of each render and concatenate them into an
%image matrix. 1
[allTopos(1:3).frames] =deal([]);%Data
prevFigH=[];
for iSubj = subjList,
    thisSubj = dirNames{iSubj};
    
    
    figH(1)= figure(300+iSubj);
    clf;
    %plotOnEgi(unstacked(:,iSubj,iT));
    plotContourOnScalp(unstackedData(:,iSubj,iT),thisSubj,projectDir)
    figH(2)=figure(400+iSubj);
    clf;
    plotContourOnScalp(unstackedYhatMN(:,iSubj,iT),thisSubj,projectDir);
     
    figH(3)=figure(500+iSubj);
    clf;
    plotContourOnScalp(unstackedYhatLASSO(:,iSubj,iT),thisSubj,projectDir);
     
    figTitles = {...
    {'Data';['Participant: ' thisSubj]}...
    {'Minimum-Norm Fit';['Participant: ' thisSubj]}...
    {'Group LASSO fit';['Participant: ' thisSubj]} };
    

    for iFig = 1:3,
        figure(figH(iFig))
        view(20,35)
        camproj('perspective')
        axis off;
        set(gcf,'position',pos+(iFig-1)*[362 0 0 0]);
        set(gcf,'color',[1 1 1]);
        title(figTitles{iFig});
        thisFrame = getframe(figH(iFig),[90 100 190 300]);
        allTopos(iFig).frames = [allTopos(iFig).frames thisFrame.cdata];
    end
    
    %if rendering multiple topos close ones we don't need. 
    if closeTopoFigs
        if ~isempty(prevFigH)
            close(prevFigH)
        end
        prevFigH = figH;
        
    end
    drawnow;
  
end

%%
%Plot all topoggraphies on a single figure. 
if strcmp(numChoice,'All')
    figure(600)
    allTopoCat = cat(1,allTopos(:).frames);
    imagesc(allTopoCat)
    axis ij
    axis tight
    axis off
    fullScreen = get(0,'ScreenSize');
    set(gcf,'Position',fullScreen,'color','w')
 
    axis equal
    
    xPos = size(allTopos(1).frames,1);
    text(-75,round(xPos/2),'Data','fontsize',18)
    text(-75,round(xPos/2)+xPos,{'Minimum';'Norm'},'fontsize',18)
    text(-75,round(xPos/2)+2*xPos,{'Group';'Lasso'},'fontsize',18)
end


% Don't plot anything if not asked. 
    case 'No',
        
end

%% Create data topographies from each participant.

%The dimensionality reduced data is fit in the main algorithm.  Uncomment
%the code below to show topographies for the full dimension data.

%Time index to plot. 196th sample = 251.3 ms.
% iT = 196;
% 
% pos =  [80   300   362   490]
% unstacked = reshape(Y,128,9,780);
% nSubj = size(unstacked,2);
% for iSubj = 1:nSubj,
%     figure(200+iSubj);
%     clf
%     %plotOnEgi(unstacked(:,iSubj,iT));
%     plotContourOnScalp(unstacked(:,iSubj,iT),'skeri0044',projectDir)
%     view(20,35)
%     camproj('perspective')
%     axis off
%     set(gcf,'position',pos);
%     name = ['prettyScalp_subj_' num2str(iSubj) '_time_' num2str(round(iT*Axx.dTms))];
% end

%% Figure showing time courses from minimum-norm solution for every ROI on single plot.

figH = [];
leftIdx = cell2mat(arrayfun(@(x) ~isempty(strfind(ROIs{1}.name{x},'-L')), 1:length(ROIs{1}.name),'uni',false));
fontSize = 12;
gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Arial',...
    'linewidth',1,'YLim',[-2 2]*1e-4};

tempColors = colormap(jet);
tempIdx = round(linspace(0,length(tempColors),length(find(leftIdx==1))+1));
roiColors = tempColors(tempIdx(2:end),:);

t = 0:Axx.dTms:(Axx.nT-1)*Axx.dTms;

figH(1) = figure;
subplot(2,1,1);
hold on;
selectedRois = find(leftIdx==1);

selectedColors = roiColors(ceil((selectedRois)/2),:);

plot(t,regionActivityMinNorm(selectedRois,:)','linewidth',2)

xlim([0,size(regionActivityMinNorm,2)]);
legend(ROIs{1}.name(leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;
title('Min Norm Solution - Left Hemisphere');
xlabel('Time (ms)')
ylabel('Current Source Density')

subplot(2,1,2);
hold on;
selectedRois = find(leftIdx==0); %Select left ROI responses

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(t,regionActivityMinNorm(selectedRois,:)','linewidth',2)


xlim([0,size(regionActivityMinNorm,2)]);
legend(ROIs{1}.name(~leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;

title('Min Norm Solution - Right Hemisphere');
xlabel('Time (ms)')
ylabel('Current Source Density')
%
figH(2) = figure;
subplot(2,1,1);
hold on;
selectedRois = find(leftIdx==1);

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(t,regionActivity(selectedRois,:)','linewidth',2)

legend(ROIs{1}.name(leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off;

title('Group Lasso Solution - Left Hemisphere');
xlabel('Time (ms)')
ylabel('Current Source Density')


subplot(2,1,2);
hold on;
selectedRois = find(leftIdx==0);

selectedColors = roiColors(ceil((selectedRois)/2),:);
set(gca, 'ColorOrder', selectedColors, 'NextPlot', 'replacechildren');

plot(t,regionActivity(selectedRois,:)','linewidth',2)

legend(ROIs{1}.name(~leftIdx),'location','southeastoutside');
set(gca,gcaOpts{:})
hold off; 
pos = get(figH(1), 'Position');
pos(3) = pos(3)*2; % Select the height of the figure in [cm]
set(figH, 'Position', pos);
title('Group Lasso Solution - Right Hemisphere');
xlabel('Time (ms)')
ylabel('Current Source Density')

