%mse plot
gca = subplot(3,2,1);
semilogy(rsquaredRidge, mseRidge, '.', rsquared, mse, 'x', rsquaredOls, mseOls, 'o');
hold, plot([rsquaredMinNorm, rsquaredMinNorm], ylim, 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], ylim, 'green');
title(['mse - ', num2str(numSubjects), ' subject']);

%auc far plot
gca = subplot(3,2,2);
plot(rsquaredRidge, aucFarRidge, '.', rsquared, aucFar, 'x', rsquaredOls, aucFarOls, 'o');
set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1);
hold, plot([rsquaredMinNorm, rsquaredMinNorm], [0, 1], 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], [0, 1], 'green');
title('auc far');

%auc close plot
gca = subplot(3,2,3);
plot(rsquaredRidge, aucCloseRidge, '.', rsquared, aucClose, 'x', rsquaredOls, aucCloseOls, 'o');
set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1);
hold, plot([rsquaredMinNorm, rsquaredMinNorm], [0, 1], 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], [0, 1], 'green');
title('auc close');

%average both aucs
gca = subplot(3,2,4);
plot(rsquaredRidge, (aucCloseRidge+aucFarRidge)/2, '.', rsquared, (aucFar+aucClose)/2, 'x', rsquaredOls, (aucFarOls+aucCloseOls)/2, 'o');
set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1);
hold, plot([rsquaredMinNorm, rsquaredMinNorm], [0, 1], 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], [0, 1], 'green');
title('average of auc close and far');

%plot relative energy
gca = subplot(3,2,5);
plot(rsquaredRidge, energyRidge, '.', rsquared, energy, 'x', rsquaredOls, energyOls, 'o');
hold, plot([rsquaredMinNorm, rsquaredMinNorm], ylim, 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], ylim, 'green');
legend('ridge', 'lasso');
title('relative energy');

%gcv error plot
subplot(3,2,6);
gca = plot(rsquaredRidge, gcvErrorRidge, rsquared, gcvError, rsquaredOls, gcvErrorOls, 'o');
legend('ridge', 'lasso', 'Location', 'NorthEast');
hold, plot([rsquaredMinNorm, rsquaredMinNorm], ylim, 'blue');
plot([rsquared(bestIndex), rsquared(bestIndex)], ylim, 'green');
title('gcv error');
xlabel('rsquared');