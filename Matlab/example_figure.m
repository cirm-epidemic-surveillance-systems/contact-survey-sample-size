clear
close all


% Model parameters



% Assortativity parameter array of values and indices of values to plot
Alpha_arr = 0:2:40;

% Array of sigma values to use
Sigma_arr = 0:0.05:1;

% Number of bins to discretise and plot the contact matrix
nBins = 100;

% Convergence tolerance and relaxation factor for Tom's iterative method
TOL = 1e-10;
relFact = 0.5;

% Make a vector containing the proportion of the population in each
% activity level bin
pPop = (1/nBins)*ones(1, nBins);


% Grid of activity level quantiles (at bin midpoints)
c = [0, cumsum(pPop)];
x = 0.5*(c(1:end-1)+c(2:end));




% Vary alpha and sigma
na = length(Sigma_arr);
nb = length(Alpha_arr);
domEig = zeros(na, nb);

for ia = 1:na
    Sigma = Sigma_arr(ia);
    if Sigma > 0
        v = logninv(x, 0, Sigma);
    else
        % If sigma=0, everyone's activity level is the same
        v = ones(size(x));
    end

    for ib = 1:nb
        Alpha = Alpha_arr(ib);
            
        % Use Tom's method to make contract matrix according to assortativity model
        M = makeContactMatrix_AM_Tom(v, pPop, Alpha, TOL, relFact);

        % Calculate dominant eigenvalue
        domEig(ia, ib) = eigs(M, 1);
    end
end



%  Heat map of the dominant eigenvalue against epsilon and b
h = figure(1);
h.Position = [ 680   482   560   496];
imagesc(Sigma_arr, Alpha_arr, domEig'); c = colorbar
ha = gca; ha.YDir = 'normal';
xline(0, 'Color', 0.7*[1 1 1], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'no heterogeneity')
xline(0.45, 'r--', 'LineWidth', 2, 'DisplayName', 'France data')
legend('Location', 'southoutside');
xlabel('\sigma')
ylabel('\alpha')
title('dominant eigenvalue (relative to \sigma=0 case)')   

