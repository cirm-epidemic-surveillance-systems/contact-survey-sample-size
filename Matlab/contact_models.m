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
b_arr = 0:2:40;
b_toPlot = [1, 11, 21];

% Array of sigma values to use
Sigma_arr = 0:0.1:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of bins to discretise and plot the contact matrix
nBins = 100;

% Number of different values of epsilon and b
na = length(eps_arr);
nb = length(b_arr);

TOL = 1e-10;
relFact = 0.5;

% Create a grid of nBins x values placed at the midpoint of each "bin"
dx = 1/nBins;
x = dx/2:dx:(1-dx/2);
nx = length(x);

% Calculate v as the inverse CDF of the log normal activity level distribution
v = logninv(x, 0, Sigma);

% Calculate approximate moments of the distribution
Ev = dx*sum(v);

% Define matrices of x and y values for calculating M(x,y)
[X, Y] = meshgrid(x, x);


% Proportionate mixing model matrix
M_PM = dx * v'.*v/Ev;

% Calculate dominant eigenvalue for proportionate mixing
domEig_PM = eigs(M_PM, 1);

% Calculate the aggregate contacts for someone as a function of their
% activity class by summing columns of the matrix
aggCont_PM = sum(M_PM, 1);
       
   

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
        b = b_arr(ib);
            
        % Make contract matrix according to assortativity model
        M1 = makeContactMatrix(X, Y, v, b );

        % Linear combination of proportionate and assortative matrices
        M_AM = (1-eps)*M_PM + eps * M1;

        % Calculate dominant eigenvalue
        domEig(ia, ib) = eigs(M_AM, 1);

        % Calculate the aggregate contacts for someone as a function of their
        % activity class by summing columns of the matrix
        aggCont_AM = sum(M_AM, 1);


        % Tom's method
        gk = calcKernel(Y-X, b);
        w = ones(size(v));
        convFlag = false;
        while ~convFlag
            wSav = w;
            w = (1-relFact)*w + relFact * 1/dx * v./(gk*w')';
            convFlag = norm(w-wSav, inf)/norm(wSav, inf) < TOL;
        end
        T = dx* w.*w'.*gk;

        % Check matrix is symmetric to within tolerance
        assert(max(max(abs(T-T'))) < 1e-12);

        T = (1-eps)*M_PM + eps*T;
        
        domEig_Tom(ia, ib) = eigs(T, 1);
        aggCont_Tom = sum(T, 1);



        % Make plots (for selected parameter values)
        if ismember(ia, eps_toPlot) & ismember(ib, b_toPlot)
            % Plot contact matrix
            figure(1);
            nexttile;
            imagesc(x, x, M_AM);
            title(sprintf('eps = %.1f, b = %.1f, lambda = %.2f', eps, b, domEig(ia, ib) ))         
            h = gca;
            h.YDir = 'normal';
            colorbar;
            xlabel('activity level quantile of individual (x)')
            ylabel('activity level quantile of contact (y)')

            figure(2);
            nexttile;
            imagesc(x, x, T);
            title(sprintf('eps = %.1f, b = %.1f, lambda = %.2f', eps, b, domEig_Tom(ia, ib) ))         
            h = gca;
            h.YDir = 'normal';
            colorbar;
            xlabel('activity level quantile of individual (x)')
            ylabel('activity level quantile of contact (y)')


            % Plot activity level distribution 
            figure(3);
            nexttile;
            plot(x, aggCont_PM, x, aggCont_AM, x, aggCont_Tom)
            hold on 
            plot(x, v, '--')
            grid on
            xlabel('activity level quantile')
            ylabel('total contact rate')
            legend('PM', 'AM', 'Tom', 'log-norm dist', 'location', 'northwest')
            title(sprintf('eps = %.1f, b = %.1f', eps, b))
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
        Ev = dx*sum(v);
    else
        % If sigma=0, everyone's activity level is 1
        v = ones(size(x));
        Ev = 1;
    end

    % Proportionate mixing model matrix
    M_PM = v'.*v/Ev;

    for ib = 1:nb
        b = b_arr(ib);
            
        % Make contract matrix according to assortativity model
        M_AM = makeContactMatrix(X, Y, v, b );

        % Calculate dominant eigenvalue
        domEig2(ia, ib) = eigs(M_AM, 1);

        % Tom's method
        gk = calcKernel(Y-X, b);
        w = ones(size(v));
        convFlag = false;
        while ~convFlag
            wSav = w;
            w = (1-relFact)*w + relFact * 1/dx * v./(gk*w')';
            convFlag = norm(w-wSav, inf)/norm(wSav, inf) < TOL;
        end
        T = dx * w.*w'.*gk;


        % Calculate dominant eigenvalue
        domEig2_Tom(ia, ib) = eigs(T, 1);
    end
end



% Contour plot of the dominant eigenvalue against epsilon and b
h = figure(4);
h.Position = [  31         105        1194         800];
tiledlayout(2, 2, "TileSpacing", "compact")
nexttile;
contourf(b_arr, eps_arr, domEig);
h = gca;
h.YDir = 'normal';
colorbar;
clim([1.25 1.75])
xlabel('b')
ylabel('\epsilon')
title('Mike - dominant eigenvalue (relative to \sigma=0)')



% Plot of the dominant eigenvalue against sigma for various b (every 5th
% value)
bPick = 1:5:nb;
nexttile;
plot(Sigma_arr, domEig2(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue')
l = legend(string(b_arr(bPick)), 'Location', 'northwest');
title(l, 'b');

nexttile;
contourf(b_arr, eps_arr, domEig_Tom);
h = gca;
h.YDir = 'normal';
colorbar;
clim([1.25 1.75])
xlabel('b')
ylabel('\epsilon')
title('Tom - dominant eigenvalue (relative to \sigma=0)')



% Plot of the dominant eigenvalue against sigma for various b (every 5th
% value)
bPick = 1:5:nb;
nexttile;
plot(Sigma_arr, domEig2_Tom(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue')
l = legend(string(b_arr(bPick)), 'Location', 'northwest');
title(l, 'b');
