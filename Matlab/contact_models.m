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
M_PM = v'.*v/Ev;

% Calculate dominant eigenvalue for proportionate mixing
domEig_PM = dx*eigs(M_PM, 1);

% Calculate the aggregate contacts for someone as a function of their
% activity class by summing columns of the matrix
aggCont_PM = dx*sum(M_PM, 1);
       
   

% Initialise array for dominant eigenvalue
domEig = zeros(na, nb);


% Set up figure for plotting matrices
h = figure(1);
h.Position = [      93         113        1096         837];
tiledlayout(3, 3, "TileSpacing", "compact")

h = figure(2);
h.Position = [      180          94        1096         837];
tiledlayout(3, 3, "TileSpacing", "compact")

% Calculate the fully assortative mixing model (epsilon=1) for each value of b
for ia = 1:na
    eps = eps_arr(ia);
    for ib = 1:nb
        b = b_arr(ib);
            
        % Calculate the kernel as a fuction of Y-X
        gk = exp(-b*(Y-X).^2);
    
        % Calculate the product v(y)*g(y-x)
        C = v'.*gk;
    
        % Calculte denominator of equation for M
        den = dx*sum(C, 1);
    
        % Assortative mixing matrix
        M_AM = (1-eps)*M_PM + eps * v.*C./den;

        % Calculate dominant eigenvalue
        domEig(ia, ib) = dx*eigs(M_AM, 1);

        % Calculate the aggregate contacts for someone as a function of their
        % activity class by summing columns of the matrix
        aggCont_AM = dx*sum(M_AM, 1);

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

            % Plot activity level distribution 
            figure(2);
            nexttile;
            plot(x, aggCont_PM, x, aggCont_AM)
            hold on 
            plot(x, v, '--')
            grid on
            xlabel('activity level quantile')
            ylabel('total contact rate')
            legend('PM', 'AM', 'log-norm dist', 'location', 'northwest')
            title(sprintf('eps = %.1f, b = %.1f', eps, b))
        end
    end
end
sgtitle('check that matrix column sums gives correct activity level dist')


% Now do the same thing with the simplified assortative mixing model (fixing eps=1), varying b and sigma
na = length(Sigma_arr);
domEig2 = zeros(na, nb);
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
            
        % Calculate the kernel as a fuction of Y-X
        gk = exp(-b*(Y-X).^2);
    
        % Calculate the product v(y)*g(y-x)
        C = v'.*gk;
    
        % Calculte denominator of equation for M
        den = dx*sum(C, 1);
    
        % Assortative mixing matrix
        M_AM = v.*C./den;

        % Calculate dominant eigenvalue
        domEig2(ia, ib) = dx*eigs(M_AM, 1);
    end
end



% Contour plot of the dominant eigenvalue against epsilon and b
figure(3);
contourf(b_arr, eps_arr, domEig);
h = gca;
h.YDir = 'normal';
colorbar;
xlabel('b')
ylabel('\epsilon')
title('dominant eigenvalue (relative to \sigma=0)')



% Plot of the dominant eigenvalue against sigma for various b (every 5th
% value)
bPick = 1:5:nb;
figure(4);
plot(Sigma_arr, domEig2(:, bPick));
grid on
xlabel('\sigma')
ylabel('dominant eigenvalue')
l = legend(string(b_arr(bPick)), 'Location', 'northwest')
title(l, 'b');
