clear 
close all

% Pop size
N = 1e4;

% Number of "bins" (range)
nBins = [2:100, 200:100:2000];

% SD of activity distribution (log scale)
sigma = 1;

% Activity vector
activity = lognrnd(0, sigma, 1, N);

% Full contact matrix (one row/col per individual)
M = activity.*(activity')/sum(activity);

% Dominant eigenvalue of M ('true' matrix)
lambdaTrue = eigs(M, 1);

nTries = length(nBins);
lambdaBinned = zeros(1, nTries);
% Loop through different numbers of bins
for iTry = 1:nTries
    % Find quantiles of activity levels to define bin edges
    qt = 0:(1/nBins(iTry)):1;
    activityQuantiles = quantile(activity, qt);
    activityBinMean = zeros(1, nBins(iTry));
    % Find the mean activity level with each bin
    for iBin = 1:nBins(iTry)
        inBinFlag = activity >= activityQuantiles(iBin) & activity < activityQuantiles(iBin+1);
        activityBinMean(iBin) = mean(activity(inBinFlag)); 
    end
    % Form the binned contact matrix and find its dominant eigenvalue
    Mbinned = activityBinMean.*(activityBinMean')/sum(activityBinMean);
    lambdaBinned(iTry) = eigs(Mbinned, 1);
end


figure;
plot(nBins, lambdaBinned, '.-')
yline(lambdaTrue, 'r--')
xlabel('number of bins')
ylabel('dominant eigenvalue')
title(sprintf('dominant eigenvalue vs number of bins (proportionate mixing, pop size N=%i)', N))
legend('binned', '"true" value', 'Location', 'southeast')
grid on

