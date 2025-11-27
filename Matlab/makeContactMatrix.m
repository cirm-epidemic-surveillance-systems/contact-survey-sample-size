function M = makeContactMatrix(X, Y, v, b)

% Make a contact matrix using activity level distribution v and with
% assortativity kernel parmaeter b. X and Y are square matrices providing a
% grid of x values (quantiles of the activity level distribution) that
% define the activity level bins

% Number of bins
[nBins, ~] = size(X);

dx = X(1, 2)-X(1, 1);

% Calculate the kernel as a fuction of Y-X
gk = calcKernel(Y-X, b);

% Calculate the product v(y)*g(y-x)
C = v'.*gk;

% Calculte denominator of equation for M
den = dx*sum(tril(C));

% Calculate matrix v(y)*g(y-x) / int_x^1 [v(y)*g(y-x)*dy]
M1 = C./den;

% For the first column, multiply by v(x)
M1(:, 1) = M1(:, 1) * v(1);

% For each subsequent column (from the main diagonal down), multiply by [v(x) - int_0^x M(x',x)*dx']
for iCol = 2:nBins
    M1(iCol:end, iCol) = M1(iCol:end, iCol) * (v(iCol) - dx*sum(M1(iCol, 1:iCol-1)));
end

% Fill the upper triangle of the matrix by transposing the lower triangle
M = tril(M1) + tril(M1, -1)';


% Check matrix is symmetric to within tolerance
assert(max(max(abs(M-M'))) < 1e-12);

