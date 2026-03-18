clear
close all

% Visualise some contact matrices constructed from the "continuous operator
% model" M(x,y) described in the Overleaf draft, with a continuous
% distribution of activity levels and different degrees of assortativity

% Model parameters

% Std. dev. in the log activity level distribution
Sigma = 0.4;

% Assortativity parameter array of values and indices of values to plot
Alpha_arr = 0:2:40;
Alpha_toPlot = [1, 11, 21];

% Array of sigma values to use
Sigma_arr = 0:0.1:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of bins to discretise and plot the contact matrix
nBins = 100;

% Number of different values of Alpha
nb = length(Alpha_arr);

% Convergence tolerance and relaxation factor for Tom's iterative method
TOL = 1e-10;
relFact = 0.5;

% Make a vector containing the proportion of the population in each
% activity level bin
pPop = (1/nBins)*ones(1, nBins);
%pPop = linspace(10, 1, nBins); pPop = pPop/sum(pPop);
%pPop = [10*ones(1, 10), ones(1, nBins-10)]; pPop = pPop/sum(pPop);

% Grid of activity level quantiles (at bin midpoints)
c = [0, cumsum(pPop)];
x = 0.5*(c(1:end-1)+c(2:end));

% Calculate v as the inverse CDF of the log normal activity level distribution
v = logninv(x, 0, Sigma);

% Proportionate mixing model matrix
M_PM = makeContactMatrix_PM(v, pPop);

% Calculate dominant eigenvalue for proportionate mixing
domEig_PM = eigs(M_PM, 1);

% Calculate the aggregate contacts for someone as a function of their
% activity class by summing columns of the matrix
aggCont_PM = sum(M_PM, 1);
       

% Plot PM contact matrix
figure(10);
imagesc(x, x, M_PM);
title(sprintf('lambda = %.2f', domEig_PM ))         
h = gca;
h.YDir = 'normal';
colorbar;
xlabel('activity quantile of individual (x)')
ylabel('activity quantile of contact (y)')




% Initialise array for dominant eigenvalue
domEig = zeros(1, nb);

% Set up figure for plotting matrices
h = figure(1);
h.Position = [       99         633        1294         330];
tiledlayout(1, 3, "TileSpacing", "compact")

h = figure(2);
h.Position = [   148         479        1294         387];
tiledlayout(1, 3, "TileSpacing", "compact")

for ib = 1:nb
    Alpha = Alpha_arr(ib);
        
    % Make assortative contact matrix using Tom's method
    M = makeContactMatrix_AM_Tom(v, pPop, Alpha, TOL, relFact);

   
    % Calculate dominant eigenvalue
    domEig(ib) = eigs(M, 1);
    
    % Calculate the aggregate contacts for someone as a function of their
    % activity class by summing columns of the matrix
    aggCont = sum(M, 1);


    % Make plots (for selected parameter values)
    if ismember(ib, Alpha_toPlot)
        % Plot contact matrix
        figure(1);
        nexttile;
        imagesc(x, x, M);
        title("\alpha = " + sprintf('%.1f', Alpha) + ", \lambda = " + sprintf('%.2f', domEig(ib) ) )         
        h = gca;
        h.YDir = 'normal';
        colorbar;
        xlabel('activity quantile of individual (x)')
        ylabel('activity quantile of contact (y)')


        % Plot activity level distribution 
        figure(2);
        nexttile;
        plot(x, aggCont_PM, x, aggCont)
        hold on 
        [Ev, ~] = lognstat(0, Sigma);
        plot(x, v/Ev, '--')
        grid on
        xlabel('activity level quantile')
        ylabel('total contact rate')
        legend('PM', 'AM', 'target', 'location', 'northwest')
        title("\alpha = " + sprintf('%.1f', Alpha))
    end
end
sgtitle('check that matrix column sums gives correct activity level dist')






% Now varying b and sigma
na = length(Sigma_arr);
domEig2 = zeros(na, nb);

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
        domEig2(ia, ib) = eigs(M, 1);
    end
end



% Contour plot of the dominant eigenvalue against epsilon and b
h = figure(4);
% Plot of the dominant eigenvalue against sigma for various Alpha (every 5th
% value)
bPick = 1:5:nb;
plot(Sigma_arr, domEig2(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue (\lambda)')
l = legend(string(Alpha_arr(bPick)), 'Location', 'northwest');
title(l, '\alpha');
