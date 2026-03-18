function M = makeContactMatrix_AM_Tom(v, pPop, Alpha, TOL, relFact)

% Function to use Tom's method to make an assortative mixing contact matrix from a specified
% distribution of activity levels (in discrete bins) and assortativity parmaeter Alpha.
% The activity level distribution (defined by v) will be normalised so that
% the output matrix has mean activity level 1.
%
% USAGE: M = makeContactMatrix_AM(v, Alpha, TOL, relFact)
%
% INPUTS: v - 1 n x vector of monotonically increaing activity levels in
% each bin
%         pPop - 1 n x vector containing the proportion of the population in each
% activity level bin
%         Alpha - non-negative scalar parameter for the assortativity kernel
%         (larger Alpha means stronger assortativity and Alpha=0 should reduce to
%         proportionate mixing)
%         TOL - converganece tolerance for iterative method (suggested
%         value 1e-10)
%         relFact - relaxation factor between 0 and 1 for fixed-point iteration (suggested
%         value 0.5)
%
% OUTPUTS: M - n x n contact matrix whose (i,j) element is the average
% number of contacts an individual in bin j has with individuals in bin i

% Get the number of binss
nBins = length(v);

% Set up grid of activity level quantiles 
c = [0, cumsum(pPop)];
x = 0.5*(c(1:end-1)+c(2:end));

% Define matrices of x and y values for calculating M(x,y)
[X, Y] = meshgrid(x, x);

% Calculate the kernel as a fuction of Y-X
gk = calcKernel(Y-X, Alpha);

% Calculate mean activity level
Ev = sum(pPop.*v);

% Initial condition for fixed point iteration
w = ones(size(v));

% Tom's iterative method
convFlag = false;
while ~convFlag
    wSav = w;
    w = (1-relFact)*w + relFact * v./(Ev*gk*(pPop.*w)')';
    convFlag = norm(w-wSav, inf)/norm(wSav, inf) < TOL;
end
M = pPop'.*w.*w'.*gk;

% Check matrix is symmetric to within tolerance
aggCont = pPop.*M;
assert(max(max(abs(aggCont-aggCont'))) < 1e-12);

