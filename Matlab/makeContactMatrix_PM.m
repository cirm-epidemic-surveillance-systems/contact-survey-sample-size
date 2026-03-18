function M = makeContactMatrix_PM(v, pPop)

% Function to make a proportionate mixing contact matrix from a specified
% distribution of activity levels (in discrete bins)
% The activity level distribution (defined by v) will be normalised so that
% the output matrix has mean activity level 1.
%
% USAGE: M = makeContactMatrix_PM(v)
%
% INPUTS: v - 1 x n vector of monotonically increaing activity levels in
% each bin
%         pPop - 1 n x vector containing the proportion of the population in each
% activity level bin
%
% OUTPUTS: M - n x n contact matrix whose (i,j) element is the average
% number of contacts an individual in bin j has with individuals in bin i

% Get the number of bins
nBins = length(v);

% Ensure pPop is normalised
pPop = pPop/sum(pPop);

% Calculate mean activity level
Ev = sum(pPop.*v);

% Make matrix
M = pPop'.*v'.*v/Ev^2;
