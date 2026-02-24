function M = makeContactMatrix_AM_Tom(v, b, TOL, relFact)

% Function to use Tom's method to make an assortative mixing contact matrix from a specified
% distribution of activity levels (in discrete bins) and assortativity
% kernel parmaeter b.
%
% USAGE: M = makeContactMatrix_AM(v, b, TOL, relFact)
%
% INPUTS: v - 1 n x vector of monotonically increaing activity levels in
% each bin
%         b - non-negative scalar parameter for the assortativity kernel
%         (larger b means stronger assortativity and b=0 should reduce to
%         proportionate mixing)
%         TOL - converganece tolerance for iterative method
%         relFact - relaxation factor for fixed-point iteration
%
% OUTPUTS: M - n x n contact matrix whose (i,j) element is the average
% number of contacts an individual in bin j has with individuals in bin i

% Get the number of binss
nBins = length(v);

% dx is the spacing between activity level quantiles
dx = 1/nBins;

% Set up grid of activity level quantiles (at bin midpoints)
x = dx/2:dx:(1-dx/2);

% Define matrices of x and y values for calculating M(x,y)
[X, Y] = meshgrid(x, x);

% Calculate the kernel as a fuction of Y-X
gk = calcKernel(Y-X, b);

% Initial condition for fixed point iteration
w = ones(size(v));

% Tom's iterative method
convFlag = false;
while ~convFlag
    wSav = w;
    w = (1-relFact)*w + relFact * 1/dx * v./(gk*w')';
    convFlag = norm(w-wSav, inf)/norm(wSav, inf) < TOL;
end
M = dx* w.*w'.*gk;

% Check matrix is symmetric to within tolerance
assert(max(max(abs(M-M'))) < 1e-12);

