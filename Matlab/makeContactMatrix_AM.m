function M = makeContactMatrix_AM(v, pPop, Alpha)

% Function to make an assortative mixing contact matrix from a specified
% distribution of activity levels (in discrete bins) and assortativity
% kernel parmaeter Alpha.
%
% USAGE: M = makeContactMatrix_AM(v, pPop, Alpha)
%
% INPUTS: v - 1 n x vector of monotonically increaing activity levels in
% each bin
%         pPop - 1 n x vector containing the proportion of the population in each
% activity level bin
%         Alpha - non-negative scalar parameter for the assortativity kernel
%         (larger Alpha means stronger assortativity and Alpha=0 should reduce to
%         proportionate mixing)
%
% OUTPUTS: M - n x n contact matrix whose (i,j) element is the average
% number of contacts an individual in bin j has with individuals in bin i


% Get the number of bins
nBins = length(v);

% Ensure pPop is normalised
pPop = pPop/sum(pPop);

% Grid of activity level quantiles
c = [0, cumsum(pPop)];
x = 0.5*(c(1:end-1)+c(2:end));

% Define matrices of x and y values for calculating M(x,y)
[X, Y] = meshgrid(x, x);

% Calculate mean activity level
Ev = sum(pPop.*v);

% Calculate the kernel as a fuction of Y-X
gk = calcKernel(Y-X, Alpha);

% Calculate the product pi_i * v_i * g_ij
C = pPop'.*v'.*gk;

% Calculte denominator of equation for M
den = sum(tril(C));

% Calculate matrix pi_i * v_i * g_ij / sum_j^n [pi_k * v_k  * g_kj]
M1 = C./den;

% For the first column (j=1), multiply by v_j/Ev
M1(:, 1) = M1(:, 1) * v(1)/Ev;

% For each subsequent column (from the main diagonal down), multiply by [v_j - sum_1^(j-1) pi_k/pi_j * M_jk ]
for jCol = 2:nBins
    M1(jCol:end, jCol) = M1(jCol:end, jCol) * (v(jCol)/Ev - sum(pPop(1:jCol-1).*M1(jCol, 1:jCol-1))/pPop(jCol) );
end

% Fill the upper triangle of the matrix by transposing the lower triangle
% and multiplying by pi_i/pi_j
M = tril(M1) + tril(M1, -1)' .* (pPop'./pPop);


% Check detailed balance condition is satisfied to within tolerance
aggCont = pPop.*M;
assert(max(max(abs(aggCont-aggCont'))) < 1e-12);

