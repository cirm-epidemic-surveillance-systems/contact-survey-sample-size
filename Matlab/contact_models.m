clear
close all

% Visualise some contact matrices constructed from the "continuous operator
% model" M(x,y) described in the Overleaf draft, with a continuous
% distribution of activity levels and different degrees of assortativity

% Model parameters

% Std. dev. in the log activity level distribution
Sigma = 0.4;

% Assortativity constant array of values and indices of values to plot
eps_arr = 0:0.05:1;
eps_toPlot = [1, 11, 21];

% Inverse width of assortativity kernel array of values and indices of values to plot
Alpha_arr = 0:2:40;
Alpha_toPlot = [1, 11, 21];

% Array of sigma values to use
Sigma_arr = 0:0.1:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of bins to discretise and plot the contact matrix
nBins = 100;

% Number of different values of epsilon and Alpha
na = length(eps_arr);
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
xlabel('activity level quantile of individual (x)')
ylabel('activity level quantile of contact (y)')




% Initialise array for dominant eigenvalue
domEig = zeros(na, nb);

% Set up figure for plotting matrices
h = figure(1);
h.Position = [      93         113        1096         837];
tiledlayout(3, 3, "TileSpacing", "compact")

h = figure(2);
h.Position = [      180          94        1096         837];
tiledlayout(3, 3, "TileSpacing", "compact")

h = figure(3);
h.Position = [      180          94        1096         837];
tiledlayout(3, 3, "TileSpacing", "compact")

% Calculate the fully assortative mixing model (epsilon=1) for each value of b
for ia = 1:na
    eps = eps_arr(ia);
    for ib = 1:nb
        Alpha = Alpha_arr(ib);
            
        % Make contract matrix according to assortativity model
        M_AM = makeContactMatrix_AM(v, pPop, Alpha);

        % Linear combination of proportionate and assortative matrices
        M = (1-eps)*M_PM + eps * M_AM;

        domEig(ia, ib) = eigs(M, 1);

        % Calculate the aggregate contacts for someone as a function of their
        % activity class by summing columns of the matrix
        aggCont = sum(M, 1);


        % Make assortative contact matrix using Tom's method
        T_AM = makeContactMatrix_AM_Tom(v, pPop, Alpha, TOL, relFact);

        % Linear combination of proportionate and assortative matrices
        T = (1-eps)*M_PM + eps*T_AM;
        
        % Calculate dominant eigenvalue
        domEig_Tom(ia, ib) = eigs(T, 1);
        
        % Calculate the aggregate contacts for someone as a function of their
        % activity class by summing columns of the matrix
        aggCont_Tom = sum(T, 1);


        % Make plots (for selected parameter values)
        if ismember(ia, eps_toPlot) & ismember(ib, Alpha_toPlot)
            % Plot contact matrix
            figure(1);
            nexttile;
            imagesc(1:nBins, 1:nBins, M);
            title(sprintf('eps = %.1f, alpha = %.1f, lambda = %.2f', eps, Alpha, domEig(ia, ib) ))         
            h = gca;
            h.YDir = 'normal';
            colorbar;
            xlabel('activity level bin of individual (x)')
            ylabel('activity level bin of contact (y)')

            figure(2);
            nexttile;
            imagesc(1:nBins, 1:nBins, T);
            title(sprintf('eps = %.1f, alpha = %.1f, lambda = %.2f', eps, Alpha, domEig_Tom(ia, ib) ))         
            h = gca;
            h.YDir = 'normal';
            colorbar;
            xlabel('activity level bin of individual (x)')
            ylabel('activity level bin of contact (y)')


            % Plot activity level distribution 
            figure(3);
            nexttile;
            plot(x, aggCont_PM, x, aggCont, x, aggCont_Tom)
            hold on 
            [Ev, ~] = lognstat(0, Sigma);
            plot(x, v/Ev, '--')
            grid on
            xlabel('activity level quantile')
            ylabel('total contact rate')
            legend('PM', 'AM', 'AM(Tom)', 'target', 'location', 'northwest')
            title(sprintf('eps = %.1f, b = %.1f', eps, Alpha))
        end
    end
end
sgtitle('check that matrix column sums gives correct activity level dist')






% Now do the same thing with the simplified assortative mixing model (fixing eps=1), varying b and sigma
na = length(Sigma_arr);
domEig2 = zeros(na, nb);
domEig2_Tom = zeros(na, nb);

for ia = 1:na
    Sigma = Sigma_arr(ia);
    if Sigma > 0
        v = logninv(x, 0, Sigma);
    else
        % If sigma=0, everyone's activity level is 1
        v = ones(size(x));
    end

    % Proportionate mixing model matrix
    M_PM = makeContactMatrix_PM(v, pPop);

    for ib = 1:nb
        Alpha = Alpha_arr(ib);
            
        % Make contract matrix according to assortativity model
        M_AM = makeContactMatrix_AM(v, pPop, Alpha);

        % Calculate dominant eigenvalue
        domEig2(ia, ib) = eigs(M_AM, 1);

        % Use Tom's method to make contract matrix according to assortativity model
        T = makeContactMatrix_AM_Tom(v, pPop, Alpha, TOL, relFact);

        % Calculate dominant eigenvalue
        domEig2_Tom(ia, ib) = eigs(T, 1);
    end
end



% Contour plot of the dominant eigenvalue against epsilon and b
h = figure(4);
h.Position = [  31         105        1194         800];
tiledlayout(2, 2, "TileSpacing", "compact")
nexttile;
contourf(Alpha_arr, eps_arr, domEig);
h = gca;
h.YDir = 'normal';
colorbar;
clim([1.15 1.65])
xlabel('\alpha')
ylabel('\epsilon')
title('Mike - dominant eigenvalue (relative to \sigma=0)')



% Plot of the dominant eigenvalue against sigma for various Alpha (every 5th
% value)
bPick = 1:5:nb;
nexttile;
plot(Sigma_arr, domEig2(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue')
l = legend(string(Alpha_arr(bPick)), 'Location', 'northwest');
title(l, '\alpha');

nexttile;
contourf(Alpha_arr, eps_arr, domEig_Tom);
h = gca;
h.YDir = 'normal';
colorbar;
clim([1.15 1.65])
xlabel('\alpha')
ylabel('\epsilon')
title('Tom - dominant eigenvalue (relative to \sigma=0)')



% Plot of the dominant eigenvalue against sigma for various Alpha (every 5th
% value)
bPick = 1:5:nb;
nexttile;
plot(Sigma_arr, domEig2_Tom(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue')
l = legend(string(Alpha_arr(bPick)), 'Location', 'northwest');
title(l, '\alpha');
