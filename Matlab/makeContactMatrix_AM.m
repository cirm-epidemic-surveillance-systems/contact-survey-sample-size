function M = makeContactMatrix_AM(v, b)

% Function to make an assortative mixing contact matrix from a specified
% distribution of activity levels (in discrete bins) and assortativity
% kernel parmaeter b.
%
% USAGE: M = makeContactMatrix_AM(v, b)
%
% INPUTS: v - 1 n x vector of monotonically increaing activity levels in
% each bin
%         b - non-negative scalar parameter for the assortativity kernel
%         (larger b means stronger assortativity and b=0 should reduce to
%         proportionate mixing)
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

% Calculate the product pi(y)*v(y)*g(y-x)
C = v'.*gk;

% Calculte denominator of equation for M
den = sum(tril(C));

% Calculate matrix v(y)*g(y-x) / int_x^1 [v(y)*g(y-x)*dy]
M1 = C./den;

% For the first column, multiply by v(x)
M1(:, 1) = M1(:, 1) * v(1);

% For each subsequent column (from the main diagonal down), multiply by [v(x) - int_0^x M(x',x)*dx']
for jCol = 2:nBins
    M1(jCol:end, jCol) = M1(jCol:end, jCol) * (v(jCol) - sum(M1(jCol, 1:jCol-1)));
end

% Fill the upper triangle of the matrix by transposing the lower triangle
M = tril(M1) + tril(M1, -1)';


% Check matrix is symmetric to within tolerance
assert(max(max(abs(M-M'))) < 1e-12);

